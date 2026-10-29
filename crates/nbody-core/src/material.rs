//! Bulk materials: density and coefficient of restitution.
//!
//! Only two material properties enter Newtonian point-mass dynamics:
//!
//! * **density** fixes the radius of a body of given mass (and hence when two
//!   bodies touch, their mutual escape speed and their Roche limits);
//! * **coefficient of restitution** fixes how much kinetic energy a
//!   non-merging collision dissipates.
//!
//! Densities are *uncompressed* reference values. Real planets are denser
//! at depth because of self-compression; [`crate::structure::fit_two_layer`]
//! absorbs that by solving for the core fraction that reproduces a body's
//! observed bulk density.

use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Material {
    /// Metallic iron-nickel (planetary cores, M-type asteroids).
    Iron,
    /// Silicate rock (mantles, S-type asteroids).
    Rock,
    /// Water ice (icy moons, comet nuclei, Kuiper-belt objects).
    Ice,
    /// Liquid water.
    Water,
    /// Hydrogen-helium envelope at Jupiter's bulk density.
    Gas,
    /// Stellar plasma at the Sun's bulk density.
    Plasma,
    /// Graphite / carbonaceous material (C-type asteroids are porous versions).
    Carbon,
    /// Rubble pile: porous regolith aggregate, very dissipative.
    Regolith,
    /// A bouncy test material for playing with restitution.
    Rubber,
    /// User-defined material.
    Custom {
        name: String,
        /// kg m^-3
        density: f64,
        /// Normal coefficient of restitution in [0, 1].
        restitution: f64,
        /// RGB display colour.
        #[serde(default = "default_custom_colour")]
        colour: [u8; 3],
    },
}

fn default_custom_colour() -> [u8; 3] {
    [180, 180, 180]
}

impl Material {
    pub const BUILTIN: [Material; 9] = [
        Material::Iron,
        Material::Rock,
        Material::Ice,
        Material::Water,
        Material::Gas,
        Material::Plasma,
        Material::Carbon,
        Material::Regolith,
        Material::Rubber,
    ];

    pub fn custom(name: impl Into<String>, density: f64, restitution: f64) -> Material {
        Material::Custom {
            name: name.into(),
            density,
            restitution: restitution.clamp(0.0, 1.0),
            colour: default_custom_colour(),
        }
    }

    /// Bulk density in kg m^-3.
    pub fn density(&self) -> f64 {
        match self {
            Material::Iron => 7874.0,
            Material::Rock => 3000.0,
            Material::Ice => 917.0,
            Material::Water => 1000.0,
            Material::Gas => 1326.0,
            Material::Plasma => 1408.0,
            Material::Carbon => 2260.0,
            Material::Regolith => 1500.0,
            Material::Rubber => 1100.0,
            Material::Custom { density, .. } => *density,
        }
    }

    /// Normal coefficient of restitution `e` (ratio of separation to approach
    /// speed along the contact normal). `0` is perfectly inelastic.
    pub fn restitution(&self) -> f64 {
        match self {
            Material::Iron => 0.6,
            Material::Rock => 0.4,
            Material::Ice => 0.3,
            Material::Water => 0.0,
            Material::Gas => 0.0,
            Material::Plasma => 0.0,
            Material::Carbon => 0.5,
            Material::Regolith => 0.1,
            Material::Rubber => 0.9,
            Material::Custom { restitution, .. } => *restitution,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Material::Iron => "iron",
            Material::Rock => "rock",
            Material::Ice => "ice",
            Material::Water => "water",
            Material::Gas => "gas",
            Material::Plasma => "plasma",
            Material::Carbon => "carbon",
            Material::Regolith => "regolith",
            Material::Rubber => "rubber",
            Material::Custom { name, .. } => name,
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Material::Iron => "metallic Fe-Ni; planetary cores",
            Material::Rock => "silicate rock; mantles and crusts",
            Material::Ice => "water ice; icy moons and comets",
            Material::Water => "liquid water",
            Material::Gas => "H/He envelope at Jupiter's bulk density",
            Material::Plasma => "stellar plasma at the Sun's bulk density",
            Material::Carbon => "graphite / carbonaceous",
            Material::Regolith => "porous rubble pile",
            Material::Rubber => "bouncy test material",
            Material::Custom { .. } => "user-defined",
        }
    }

    pub fn colour(&self) -> [u8; 3] {
        match self {
            Material::Iron => [150, 150, 165],
            Material::Rock => [160, 115, 80],
            Material::Ice => [200, 232, 255],
            Material::Water => [60, 120, 220],
            Material::Gas => [230, 190, 120],
            Material::Plasma => [255, 215, 70],
            Material::Carbon => [70, 70, 70],
            Material::Regolith => [125, 115, 100],
            Material::Rubber => [210, 60, 60],
            Material::Custom { colour, .. } => *colour,
        }
    }

    /// Parse a built-in material name (case-insensitive).
    pub fn parse(s: &str) -> Option<Material> {
        let s = s.trim().to_ascii_lowercase();
        Material::BUILTIN.iter().find(|m| m.name() == s).cloned()
    }
}

impl fmt::Display for Material {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (rho={:.0} kg/m^3, e={:.2})", self.name(), self.density(), self.restitution())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_names_round_trip() {
        for m in Material::BUILTIN.iter() {
            assert_eq!(Material::parse(m.name()).as_ref(), Some(m));
        }
        assert!(Material::parse("unobtainium").is_none());
    }

    #[test]
    fn serde_unit_variants_are_plain_strings() {
        #[derive(Deserialize)]
        struct Wrap {
            material: Material,
        }
        let w: Wrap = toml::from_str(r#"material = "iron""#).unwrap();
        assert_eq!(w.material, Material::Iron);
        let w: Wrap =
            toml::from_str(r#"material = { custom = { name = "foam", density = 50, restitution = 0.2 } }"#).unwrap();
        assert_eq!(w.material.density(), 50.0);
        assert_eq!(w.material.name(), "foam");
    }
}
