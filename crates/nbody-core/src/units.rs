//! Physical constants (CODATA 2018 / IAU 2015 nominal values) and unit systems.
//!
//! The simulation core is unit-agnostic: it only ever sees a gravitational
//! constant `G` and numbers. Scenario files declare a [`UnitSystem`] and are
//! converted to SI on load, which keeps every internal quantity in metres,
//! kilograms and seconds. Choosing SI internally costs nothing in `f64`
//! (the dynamic range of a planetary system is ~1e-11 .. 1e40, far inside
//! the ~1e308 range of a double) and removes a whole class of unit bugs.

use serde::{Deserialize, Serialize};
use std::fmt;

/// Newtonian constant of gravitation, m^3 kg^-1 s^-2 (CODATA 2018).
pub const G_SI: f64 = 6.674_30e-11;

// ---- length -----------------------------------------------------------------
pub const METRE: f64 = 1.0;
pub const KILOMETRE: f64 = 1.0e3;
/// Astronomical unit (IAU 2012 exact definition).
pub const AU: f64 = 1.495_978_707e11;
pub const LIGHT_YEAR: f64 = 9.460_730_472_580_8e15;
pub const PARSEC: f64 = 3.085_677_581_491_367e16;
pub const R_SUN: f64 = 6.957e8;
pub const R_EARTH: f64 = 6.371_0e6;
pub const R_JUPITER: f64 = 6.991_1e7;
pub const R_MOON: f64 = 1.737_4e6;

// ---- mass -------------------------------------------------------------------
pub const KILOGRAM: f64 = 1.0;
pub const M_SUN: f64 = 1.988_47e30;
pub const M_EARTH: f64 = 5.972_4e24;
pub const M_JUPITER: f64 = 1.898_2e27;
pub const M_MOON: f64 = 7.342e22;

// ---- time -------------------------------------------------------------------
pub const SECOND: f64 = 1.0;
pub const MINUTE: f64 = 60.0;
pub const HOUR: f64 = 3600.0;
pub const DAY: f64 = 86_400.0;
/// Julian year (exactly 365.25 days), the IAU standard.
pub const YEAR: f64 = 365.25 * DAY;

/// Scale factors from a user-facing unit system to SI.
///
/// `length`, `mass` and `time` are the sizes of one user unit in metres,
/// kilograms and seconds. The gravitational constant expressed in these
/// units is `G_SI * mass * time^2 / length^3`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct UnitSystem {
    pub length: f64,
    pub mass: f64,
    pub time: f64,
}

impl Default for UnitSystem {
    fn default() -> Self {
        UnitSystem::SI
    }
}

impl UnitSystem {
    pub const SI: UnitSystem = UnitSystem { length: 1.0, mass: 1.0, time: 1.0 };
    /// AU, solar masses and years: `G ≈ 4π²`.
    pub const ASTRONOMICAL: UnitSystem = UnitSystem { length: AU, mass: M_SUN, time: YEAR };
    /// AU, solar masses and days (`G = k²` with Gauss' constant k ≈ 0.01720).
    pub const GAUSSIAN: UnitSystem = UnitSystem { length: AU, mass: M_SUN, time: DAY };
    /// Earth radii, Earth masses and hours: handy for satellites and moons.
    pub const GEOCENTRIC: UnitSystem = UnitSystem { length: R_EARTH, mass: M_EARTH, time: HOUR };

    pub fn new(length: f64, mass: f64, time: f64) -> Self {
        UnitSystem { length, mass, time }
    }

    /// `G` in these units (dimensionally `L^3 M^-1 T^-2`).
    pub fn g(&self) -> f64 {
        G_SI * self.mass * self.time * self.time / (self.length * self.length * self.length)
    }

    pub fn velocity(&self) -> f64 {
        self.length / self.time
    }
    pub fn acceleration(&self) -> f64 {
        self.length / (self.time * self.time)
    }
    pub fn energy(&self) -> f64 {
        self.mass * self.velocity() * self.velocity()
    }
    pub fn density(&self) -> f64 {
        self.mass / (self.length * self.length * self.length)
    }

    pub fn length_from_name(name: &str) -> Option<f64> {
        Some(match name.trim().to_ascii_lowercase().as_str() {
            "m" | "metre" | "meter" | "metres" | "meters" => METRE,
            "km" | "kilometre" | "kilometer" => KILOMETRE,
            "au" | "astronomical_unit" => AU,
            "ly" | "lightyear" | "light_year" => LIGHT_YEAR,
            "pc" | "parsec" => PARSEC,
            "kpc" | "kiloparsec" => 1e3 * PARSEC,
            "r_sun" | "solar_radius" | "rsun" => R_SUN,
            "r_earth" | "earth_radius" | "rearth" => R_EARTH,
            "r_jupiter" | "jupiter_radius" => R_JUPITER,
            _ => return None,
        })
    }

    pub fn mass_from_name(name: &str) -> Option<f64> {
        Some(match name.trim().to_ascii_lowercase().as_str() {
            "kg" | "kilogram" | "kilograms" => KILOGRAM,
            "g" | "gram" => 1e-3,
            "t" | "tonne" => 1e3,
            "sun" | "solar" | "m_sun" | "msun" | "solar_mass" => M_SUN,
            "earth" | "m_earth" | "mearth" | "earth_mass" => M_EARTH,
            "jupiter" | "m_jupiter" | "mjup" | "jupiter_mass" => M_JUPITER,
            "moon" | "m_moon" | "lunar" | "lunar_mass" => M_MOON,
            _ => return None,
        })
    }

    pub fn time_from_name(name: &str) -> Option<f64> {
        Some(match name.trim().to_ascii_lowercase().as_str() {
            "s" | "sec" | "second" | "seconds" => SECOND,
            "min" | "minute" | "minutes" => MINUTE,
            "h" | "hr" | "hour" | "hours" => HOUR,
            "d" | "day" | "days" => DAY,
            "yr" | "year" | "years" | "julian_year" => YEAR,
            "kyr" => 1e3 * YEAR,
            "myr" => 1e6 * YEAR,
            "gyr" => 1e9 * YEAR,
            _ => return None,
        })
    }

    /// Build a unit system from names such as `("au", "solar", "year")`.
    pub fn from_names(length: &str, mass: &str, time: &str) -> Result<Self, UnitError> {
        Ok(UnitSystem {
            length: Self::length_from_name(length).ok_or_else(|| UnitError::UnknownLength(length.to_string()))?,
            mass: Self::mass_from_name(mass).ok_or_else(|| UnitError::UnknownMass(mass.to_string()))?,
            time: Self::time_from_name(time).ok_or_else(|| UnitError::UnknownTime(time.to_string()))?,
        })
    }
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum UnitError {
    #[error("unknown length unit '{0}' (try m, km, au, ly, pc, r_sun, r_earth)")]
    UnknownLength(String),
    #[error("unknown mass unit '{0}' (try kg, sun, earth, jupiter, moon)")]
    UnknownMass(String),
    #[error("unknown time unit '{0}' (try s, min, hour, day, year)")]
    UnknownTime(String),
}

/// Human-friendly formatting of an SI time: picks seconds, hours, days or years.
pub fn format_time(seconds: f64) -> String {
    let a = seconds.abs();
    if a < 2.0 * MINUTE {
        format!("{seconds:.3} s")
    } else if a < 2.0 * HOUR {
        format!("{:.3} min", seconds / MINUTE)
    } else if a < 2.0 * DAY {
        format!("{:.3} h", seconds / HOUR)
    } else if a < 2.0 * YEAR {
        format!("{:.3} d", seconds / DAY)
    } else {
        format!("{:.4} yr", seconds / YEAR)
    }
}

/// Human-friendly formatting of an SI length.
pub fn format_length(metres: f64) -> String {
    let a = metres.abs();
    if a < 1e4 {
        format!("{metres:.3} m")
    } else if a < 0.01 * AU {
        format!("{:.3} km", metres / KILOMETRE)
    } else if a < 0.1 * PARSEC {
        format!("{:.5} au", metres / AU)
    } else {
        format!("{:.4} pc", metres / PARSEC)
    }
}

/// Human-friendly formatting of an SI mass.
pub fn format_mass(kg: f64) -> String {
    let a = kg.abs();
    if a < 1e-3 * M_MOON {
        format!("{kg:.4e} kg")
    } else if a < 0.1 * M_EARTH {
        format!("{:.4} M_moon", kg / M_MOON)
    } else if a < 0.1 * M_JUPITER {
        format!("{:.4} M_earth", kg / M_EARTH)
    } else if a < 0.05 * M_SUN {
        format!("{:.4} M_jup", kg / M_JUPITER)
    } else {
        format!("{:.4} M_sun", kg / M_SUN)
    }
}

/// Human-friendly formatting of an SI speed.
pub fn format_speed(mps: f64) -> String {
    if mps.abs() < 1e3 {
        format!("{mps:.3} m/s")
    } else {
        format!("{:.4} km/s", mps / 1e3)
    }
}

impl fmt::Display for UnitSystem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "length={} m, mass={} kg, time={} s, G={:.6e}", self.length, self.mass, self.time, self.g())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn astronomical_g_is_four_pi_squared() {
        // Kepler's third law in AU/M_sun/yr gives G = 4π² to ~0.1% (the Earth's
        // sidereal year is not exactly a Julian year and M_sun absorbs the planets).
        let g = UnitSystem::ASTRONOMICAL.g();
        let four_pi_sq = 4.0 * std::f64::consts::PI * std::f64::consts::PI;
        assert!((g / four_pi_sq - 1.0).abs() < 2e-3, "G = {g}");
    }

    #[test]
    fn gaussian_constant() {
        let k = UnitSystem::GAUSSIAN.g().sqrt();
        assert!((k - 0.017_202_098_95).abs() < 2e-6, "k = {k}");
    }

    #[test]
    fn names_resolve() {
        let u = UnitSystem::from_names("AU", "Solar", "Year").unwrap();
        assert_eq!(u, UnitSystem::ASTRONOMICAL);
        assert!(UnitSystem::from_names("furlong", "kg", "s").is_err());
    }
}
