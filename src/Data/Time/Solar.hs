-- | This module contains basic solar calculation functions. It is based on the
-- methods found at https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
--
-- Accordingly, the same caveats apply. This method is therefore only accurate
-- for dates between 1901 and 2099. The sunrise and sunset results are
-- theoretically accurate to within a minute for locations between +/- 72°
-- latitude, and within 10 minutes outside of those latitudes.  
module Data.Time.Solar
    ( Location(..)
    , solarNoon
    , solarMidnight
    , sunrise
    , sunset
    , hourAngle
    , trueSolarTime
    , solarZenithAngle
    , solarElevationAngle
    -- * Re-Export
    , ZonedTime(..)
    ) where

import Data.Fixed
import Data.Time
import Data.Time.LocalTime

data Location = Location
    { latitude :: Double
    , longitude :: Double
    } deriving (Show, Eq, Read, Ord)

degToRad :: Double -> Double
degToRad deg = pi * deg / 180
{-# INLINE degToRad #-}

radToDeg :: Double -> Double
radToDeg rad = rad * 180 / pi
{-# INLINE radToDeg #-}

startDay :: UTCTime
startDay = UTCTime (fromGregorian 1900 1 1) 0
{-# INLINE startDay #-}

-- | Offset in hours
timezoneOffset :: ZonedTime -> Double
timezoneOffset = (/ 60) . fromIntegral . timeZoneMinutes . zonedTimeZone

pastLocalMidnight :: ZonedTime -> DiffTime
pastLocalMidnight = timeOfDayToTime . localTimeOfDay . zonedTimeToLocalTime

julianDay :: ZonedTime -> Double
julianDay t =
    let ts = diffUTCTime (zonedTimeToUTC t) startDay
        num = (fromIntegral $ truncate ts) / 24 / 60 / 60 + 2
        pm = (/ 86400) . fromIntegral . truncate . pastLocalMidnight $ t
     in num + 2415018.5 + pm

julianCentury :: ZonedTime -> Double
julianCentury t = (julianDay t - 2451545) / 36525
{-# INLINE julianCentury #-}

-- | deg
geomMeanLongSun :: ZonedTime -> Double
geomMeanLongSun t =
    let cent = julianCentury t
     in (280.46646 + (cent * (36000.76983 + cent * 0.0003032))) `mod'` 360

-- | deg
geomMeanAnomSun :: ZonedTime -> Double
geomMeanAnomSun t =
    let cent = julianCentury t
     in 357.52911 + cent * (35999.05029 - 0.0001537 * cent)

eccentEarthOrbit :: ZonedTime -> Double
eccentEarthOrbit t =
    let cent = julianCentury t
     in 0.016708634 - cent * (0.000042037 + 0.0000001267 * cent)

sunEqOfCtr :: ZonedTime -> Double
sunEqOfCtr t =
    let cent = julianCentury t
        anom = degToRad $ geomMeanAnomSun t
     in sin anom * (1.914602 - cent * (0.004817 + 0.000014 * cent)) +
        sin (2 * anom) * (0.019993 - 0.000101 * cent) +
        sin (3 * anom) * 0.000289

-- | deg
sunTrueLong :: ZonedTime -> Double
sunTrueLong t = geomMeanLongSun t + sunEqOfCtr t
{-# INLINE sunTrueLong #-}

-- | deg
sunTrueAnom :: ZonedTime -> Double
sunTrueAnom t = geomMeanAnomSun t + sunEqOfCtr t
{-# INLINE sunTrueAnom #-}

-- | AUs
sunRadVector :: ZonedTime -> Double
sunRadVector t =
    let ec = eccentEarthOrbit t
     in (1.000001018 * (1 - ec * ec)) /
        (1 + ec * cos (degToRad (sunTrueAnom t)))

-- | deg
sunAppLong :: ZonedTime -> Double
sunAppLong t =
    sunTrueLong t - 0.00569 -
    0.00478 * sin (degToRad (125.04 - 1934.136 * julianCentury t))

-- | deg
meanObliqEcliptic :: ZonedTime -> Double
meanObliqEcliptic t =
    let cent = julianCentury t
     in 23 +
        (26 +
         ((21.448 - cent * (46.815 + cent * (0.00059 - cent * 0.001813)))) / 60) /
        60

-- | deg
obliqCorr :: ZonedTime -> Double
obliqCorr t =
    meanObliqEcliptic t +
    0.00256 * cos (degToRad (125.04 - 1934.136 * julianCentury t))

-- | deg
sunRtAscen :: ZonedTime -> Double
sunRtAscen t =
    meanObliqEcliptic t +
    0.00256 * cos (degToRad (125.04 - 1934.136 * julianCentury t))

-- | deg
sunDeclin :: ZonedTime -> Double
sunDeclin t =
    radToDeg . asin $
    sin (degToRad $ obliqCorr t) * sin (degToRad $ sunAppLong t)

eqOfTime :: ZonedTime -> Double
eqOfTime t =
    let u2 = tan (degToRad (obliqCorr t / 2)) ** 2
        i2 = degToRad $ geomMeanLongSun t
        j2 = degToRad $ geomMeanAnomSun t
        k2 = eccentEarthOrbit t
     in 4 *
        radToDeg
            (u2 * sin (2 * i2) - 2 * k2 * sin j2 +
             4 * k2 * u2 * sin j2 * cos (2 * i2) -
             0.5 * u2 * u2 * sin (4 * i2) -
             1.25 * k2 * k2 * sin (2 * j2))

-- | deg
haSunrise :: ZonedTime -> Location -> Double
haSunrise t loc =
    let x = cos . degToRad $ 90.833
        lat = degToRad . latitude $ loc
        decl = degToRad . sunDeclin $ t
     in radToDeg . acos $ x / (cos lat * cos decl) - tan lat * tan decl

solarNoonLST :: ZonedTime -> Location -> Double
solarNoonLST t loc =
    (720 - 4 * longitude loc - eqOfTime t + timezoneOffset t * 60) / 1440

sunriseLST :: ZonedTime -> Location -> Double
sunriseLST t loc = solarNoonLST t loc - haSunrise t loc * 4 / 1440

sunsetLST :: ZonedTime -> Location -> Double
sunsetLST t loc = solarNoonLST t loc + haSunrise t loc * 4 / 1440

sunlightDuration :: ZonedTime -> Location -> DiffTime
sunlightDuration t loc = fromInteger . truncate $ 60 * 8 * haSunrise t loc

lstToLocal :: LocalTime -> Double -> LocalTime
lstToLocal t x = utcToLocalTime utc . addUTCTime x' . localTimeToUTC utc $ t'
  where
    t' = LocalTime (localDay t) midnight
    x' = fromIntegral . truncate . (* 86400) $ x

-- | Return the time of solar noon in a given timezone for a given day, at a
-- given location. Solar noon is the moment when the sun contacts the
-- observer's meridian, reaching its highest position above the horizon on that
-- day.
solarNoon :: ZonedTime -> Location -> ZonedTime
solarNoon t loc =
    ZonedTime
        (lstToLocal (zonedTimeToLocalTime t) (solarNoonLST t loc))
        (zonedTimeZone t)

-- | Return the time of solar midnight, opposite of 'solarNoon'. Note that this
-- will return the /next/ solar midnight!
solarMidnight :: ZonedTime -> Location -> ZonedTime
solarMidnight t loc =
    ZonedTime
        (lstToLocal (zonedTimeToLocalTime t) (solarNoonLST t loc + 0.5))
        (zonedTimeZone t)

-- | Determine the time of sunrise relative to a zoned time, at a given location.
sunrise :: ZonedTime -> Location -> ZonedTime
sunrise t loc =
    ZonedTime
        (lstToLocal (zonedTimeToLocalTime t) (sunriseLST t loc))
        (zonedTimeZone t)

-- | Determine the time of sunset relative to a zoned time, at a given location.
sunset :: ZonedTime -> Location -> ZonedTime
sunset t loc =
    ZonedTime
        (lstToLocal (zonedTimeToLocalTime t) (sunsetLST t loc))
        (zonedTimeZone t)

trueSolarTime' :: ZonedTime -> Location -> Double
trueSolarTime' t loc =
    let pm = (/ 60) . fromIntegral . truncate . pastLocalMidnight $ t
     in (pm + eqOfTime t + 4 * longitude loc - 60 * latitude loc) `mod'` 1440

-- | Time of day at a given location as measured by the movement of the sun,
-- given as time after midnight.
trueSolarTime :: ZonedTime -> Location -> DiffTime
trueSolarTime t loc = fromInteger . truncate $ trueSolarTime' t loc

-- | Return hour angle, one of the coordinates used in the equatorial
-- coordinate system to give the direction of a point on the celestial sphere.
-- Given in /degrees/.
hourAngle :: ZonedTime -> Location -> Double
hourAngle t loc
    | tst < 0 = tst + 180
    | otherwise = tst - 180
  where
    tst = trueSolarTime' t loc / 4

-- | The solar zenith angle is the angle between the zenith and the centre of
-- the Sun's disc. Given in /degrees/.
solarZenithAngle :: ZonedTime -> Location -> Double
solarZenithAngle t loc =
    let sins = sin (degToRad . latitude $ loc) * sin (degToRad . sunDeclin $ t)
        coss =
            cos (degToRad . latitude $ loc) * cos (degToRad . sunDeclin $ t) *
            cos (degToRad $ hourAngle t loc)
     in radToDeg . acos $ sins + coss

-- | The solar elevation angle is the altitude of the Sun, the angle between
-- the horizon and the centre of the Sun's disc. Given in /degrees/.
-- Complimentary to 'solarZenithAngle'.
solarElevationAngle :: ZonedTime -> Location -> Double
solarElevationAngle t loc = 90 - solarZenithAngle t loc
{-# INLINE solarElevationAngle #-}
