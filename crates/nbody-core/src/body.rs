//! A gravitating body: name, kinematic state and layered structure.

use crate::material::Material;
use crate::structure::{SphericalStructure, StructureError};
use crate::vector::Vec3;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Stable identity of a body inside one [`crate::system::NBodySystem`].
/// Ids survive removals and mergers (the survivor keeps its id), which is
/// what lets a projectile be tracked through an encounter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct BodyId(pub u64);

impl BodyId {
    pub const UNASSIGNED: BodyId = BodyId(0);
}

impl fmt::Display for BodyId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "#{}", self.0)
    }
}

/// A purely descriptive label; it never changes the physics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum BodyKind {
    Star,
    #[default]
    Planet,
    Moon,
    Asteroid,
    Comet,
    Projectile,
    /// Product of a merger.
    Remnant,
    /// Test-mass-like particle in a cluster.
    Particle,
}

impl BodyKind {
    pub fn name(self) -> &'static str {
        match self {
            BodyKind::Star => "star",
            BodyKind::Planet => "planet",
            BodyKind::Moon => "moon",
            BodyKind::Asteroid => "asteroid",
            BodyKind::Comet => "comet",
            BodyKind::Projectile => "projectile",
            BodyKind::Remnant => "remnant",
            BodyKind::Particle => "particle",
        }
    }

    pub fn parse(s: &str) -> Option<BodyKind> {
        Some(match s.trim().to_ascii_lowercase().as_str() {
            "star" => BodyKind::Star,
            "planet" => BodyKind::Planet,
            "moon" => BodyKind::Moon,
            "asteroid" => BodyKind::Asteroid,
            "comet" => BodyKind::Comet,
            "projectile" => BodyKind::Projectile,
            "remnant" => BodyKind::Remnant,
            "particle" => BodyKind::Particle,
            _ => return None,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Body {
    id: BodyId,
    pub name: String,
    pub kind: BodyKind,
    /// Metres.
    pub position: Vec3,
    /// Metres per second.
    pub velocity: Vec3,
    structure: SphericalStructure,
    /// Cached from `structure`.
    mass: f64,
    /// Cached from `structure`.
    radius: f64,
}

impl Body {
    pub fn new(name: impl Into<String>, structure: SphericalStructure, position: Vec3, velocity: Vec3) -> Body {
        let mass = structure.mass();
        let radius = structure.radius();
        Body { id: BodyId::UNASSIGNED, name: name.into(), kind: BodyKind::default(), position, velocity, structure, mass, radius }
    }

    /// A homogeneous body of the given material and mass.
    pub fn uniform(
        name: impl Into<String>,
        material: Material,
        mass: f64,
        position: Vec3,
        velocity: Vec3,
    ) -> Result<Body, StructureError> {
        Ok(Body::new(name, SphericalStructure::uniform_with_mass(material, mass)?, position, velocity))
    }

    pub fn with_kind(mut self, kind: BodyKind) -> Body {
        self.kind = kind;
        self
    }

    pub fn id(&self) -> BodyId {
        self.id
    }

    pub(crate) fn set_id(&mut self, id: BodyId) {
        self.id = id;
    }

    #[inline]
    pub fn mass(&self) -> f64 {
        self.mass
    }

    #[inline]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    pub fn structure(&self) -> &SphericalStructure {
        &self.structure
    }

    pub fn set_structure(&mut self, structure: SphericalStructure) {
        self.mass = structure.mass();
        self.radius = structure.radius();
        self.structure = structure;
    }

    /// Change the mass by homologous rescaling (densities fixed, radius ∝ M^{1/3}).
    pub fn set_mass(&mut self, mass: f64) -> Result<(), StructureError> {
        let s = self.structure.scaled_to_mass(mass)?;
        self.set_structure(s);
        Ok(())
    }

    /// Change the radius by homologous rescaling (densities fixed, mass ∝ R³).
    pub fn set_radius(&mut self, radius: f64) -> Result<(), StructureError> {
        let s = self.structure.scaled_to_radius(radius)?;
        self.set_structure(s);
        Ok(())
    }

    /// Replace the structure by a uniform sphere of `material` at fixed mass.
    pub fn set_material(&mut self, material: Material) -> Result<(), StructureError> {
        let s = self.structure.with_material(material)?;
        self.set_structure(s);
        Ok(())
    }

    /// The material carrying most of the mass.
    pub fn material(&self) -> &Material {
        self.structure.dominant_material()
    }

    pub fn mean_density(&self) -> f64 {
        self.structure.mean_density()
    }

    pub fn restitution(&self) -> f64 {
        self.structure.restitution()
    }

    pub fn colour(&self) -> [u8; 3] {
        self.structure.surface_material().colour()
    }

    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.velocity.norm_sq()
    }

    pub fn momentum(&self) -> Vec3 {
        self.velocity * self.mass
    }

    /// Orbital angular momentum about `origin`.
    pub fn angular_momentum(&self, origin: Vec3) -> Vec3 {
        (self.position - origin).cross(self.velocity) * self.mass
    }

    /// Surface escape speed `sqrt(2GM/R)`.
    pub fn escape_velocity(&self, g: f64) -> f64 {
        (2.0 * g * self.mass / self.radius).sqrt()
    }

    /// Surface gravity `GM/R²`.
    pub fn surface_gravity(&self, g: f64) -> f64 {
        g * self.mass / (self.radius * self.radius)
    }
}

impl fmt::Display for Body {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} [{}] m={:.4e} kg R={:.4e} m {} r={} v={}",
            self.name,
            self.kind.name(),
            self.mass,
            self.radius,
            self.material().name(),
            self.position,
            self.velocity
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn set_mass_rescales_radius() {
        let mut b = Body::uniform("x", Material::Rock, 1e15, Vec3::ZERO, Vec3::ZERO).unwrap();
        let r0 = b.radius();
        b.set_mass(8e15).unwrap();
        assert_relative_eq!(b.radius(), 2.0 * r0, max_relative = 1e-12);
        assert_relative_eq!(b.mass(), 8e15, max_relative = 1e-12);
    }

    #[test]
    fn set_material_keeps_mass() {
        let mut b = Body::uniform("x", Material::Rock, 1e15, Vec3::ZERO, Vec3::ZERO).unwrap();
        b.set_material(Material::Ice).unwrap();
        assert_relative_eq!(b.mass(), 1e15, max_relative = 1e-12);
        assert_eq!(b.material(), &Material::Ice);
        assert!(b.radius() > (1e15 / 3000.0 * 3.0 / (4.0 * std::f64::consts::PI)).cbrt());
    }
}
