//! Spherically symmetric, layered internal structure.
//!
//! A body is modelled as concentric shells, each of uniform [`Material`].
//! By Newton's shell theorem the external gravitational field of any
//! spherically symmetric mass distribution is identical to that of a point
//! mass, so the layering never changes the orbital dynamics. It does fix
//!
//! * the **radius** for a given mass (collision cross-section, mutual
//!   escape speed, Roche limits);
//! * the **moment of inertia** `I = (8π/15) Σ ρ_k (r_k⁵ − r_{k−1}⁵)`,
//!   reported as the dimensionless factor `I/(MR²)` (0.4 for a uniform
//!   sphere, 0.33 for Earth, 0.07 for the Sun);
//! * the **interior gravity profile** via the enclosed mass;
//! * what happens on **merger**: material is re-stacked densest-first,
//!   the end state of gravitational differentiation.

use crate::material::Material;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Layer {
    /// Outer radius of this shell in metres.
    pub outer_radius: f64,
    pub material: Material,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "Vec<Layer>", into = "Vec<Layer>")]
pub struct SphericalStructure {
    /// Ordered inside-out; radii strictly increasing.
    layers: Vec<Layer>,
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum StructureError {
    #[error("a body needs at least one layer")]
    NoLayers,
    #[error("layer radii must be positive, finite and strictly increasing outward (layer {index} has outer radius {radius})")]
    BadRadius { index: usize, radius: f64 },
    #[error("material density must be positive and finite (got {0})")]
    BadDensity(f64),
    #[error("mass must be positive and finite (got {0})")]
    BadMass(f64),
}

#[inline]
pub fn sphere_volume(r: f64) -> f64 {
    4.0 / 3.0 * PI * r * r * r
}

#[inline]
pub fn radius_of_volume(v: f64) -> f64 {
    (3.0 * v / (4.0 * PI)).cbrt()
}

impl TryFrom<Vec<Layer>> for SphericalStructure {
    type Error = StructureError;
    fn try_from(layers: Vec<Layer>) -> Result<Self, StructureError> {
        SphericalStructure::from_layers(layers)
    }
}

impl From<SphericalStructure> for Vec<Layer> {
    fn from(s: SphericalStructure) -> Vec<Layer> {
        s.layers
    }
}

impl SphericalStructure {
    /// A homogeneous sphere of the given radius.
    pub fn uniform_with_radius(material: Material, radius: f64) -> Result<Self, StructureError> {
        Self::from_layers(vec![Layer { outer_radius: radius, material }])
    }

    /// A homogeneous sphere whose radius follows from mass and density.
    pub fn uniform_with_mass(material: Material, mass: f64) -> Result<Self, StructureError> {
        if !(mass > 0.0) || !mass.is_finite() {
            return Err(StructureError::BadMass(mass));
        }
        let rho = material.density();
        if !(rho > 0.0) || !rho.is_finite() {
            return Err(StructureError::BadDensity(rho));
        }
        Self::uniform_with_radius(material, radius_of_volume(mass / rho))
    }

    /// Build from shells; they are sorted by radius and validated.
    pub fn from_layers(mut layers: Vec<Layer>) -> Result<Self, StructureError> {
        if layers.is_empty() {
            return Err(StructureError::NoLayers);
        }
        layers.sort_by(|a, b| a.outer_radius.partial_cmp(&b.outer_radius).unwrap_or(std::cmp::Ordering::Equal));
        let mut prev = 0.0;
        for (index, l) in layers.iter().enumerate() {
            if !(l.outer_radius > prev) || !l.outer_radius.is_finite() {
                return Err(StructureError::BadRadius { index, radius: l.outer_radius });
            }
            let rho = l.material.density();
            if !(rho > 0.0) || !rho.is_finite() {
                return Err(StructureError::BadDensity(rho));
            }
            prev = l.outer_radius;
        }
        Ok(SphericalStructure { layers })
    }

    pub fn layers(&self) -> &[Layer] {
        &self.layers
    }

    /// Outer radius in metres.
    pub fn radius(&self) -> f64 {
        self.layers.last().map(|l| l.outer_radius).unwrap_or(0.0)
    }

    pub fn volume(&self) -> f64 {
        sphere_volume(self.radius())
    }

    /// Inner radius of layer `k` (0 for the core).
    fn inner_radius(&self, k: usize) -> f64 {
        if k == 0 {
            0.0
        } else {
            self.layers[k - 1].outer_radius
        }
    }

    pub fn mass(&self) -> f64 {
        self.layers
            .iter()
            .enumerate()
            .map(|(k, l)| l.material.density() * (sphere_volume(l.outer_radius) - sphere_volume(self.inner_radius(k))))
            .sum()
    }

    pub fn mean_density(&self) -> f64 {
        self.mass() / self.volume()
    }

    /// Moment of inertia about any axis through the centre, kg m².
    pub fn moment_of_inertia(&self) -> f64 {
        self.layers
            .iter()
            .enumerate()
            .map(|(k, l)| {
                let ro = l.outer_radius;
                let ri = self.inner_radius(k);
                8.0 * PI / 15.0 * l.material.density() * (ro.powi(5) - ri.powi(5))
            })
            .sum()
    }

    /// Normalised moment of inertia `I / (M R²)`; 0.4 for a uniform sphere,
    /// smaller for centrally condensed bodies.
    pub fn inertia_factor(&self) -> f64 {
        let r = self.radius();
        self.moment_of_inertia() / (self.mass() * r * r)
    }

    /// Mass enclosed within radius `r` (clamped to the total mass outside).
    pub fn enclosed_mass(&self, r: f64) -> f64 {
        let mut m = 0.0;
        for (k, l) in self.layers.iter().enumerate() {
            let ri = self.inner_radius(k);
            if r <= ri {
                break;
            }
            let ro = l.outer_radius.min(r);
            m += l.material.density() * (sphere_volume(ro) - sphere_volume(ri));
        }
        m
    }

    /// Gravitational acceleration at radius `r` (inside or outside), m s⁻².
    pub fn gravity_at(&self, r: f64, g: f64) -> f64 {
        if r <= 0.0 {
            0.0
        } else {
            g * self.enclosed_mass(r) / (r * r)
        }
    }

    pub fn surface_material(&self) -> &Material {
        &self.layers.last().expect("validated non-empty").material
    }

    /// The material carrying the most mass.
    pub fn dominant_material(&self) -> &Material {
        let comp = self.composition();
        let best = comp
            .iter()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(m, _)| m.clone())
            .expect("validated non-empty");
        self.layers.iter().map(|l| &l.material).find(|m| **m == best).expect("present")
    }

    /// Mass per distinct material, inside-out order of first appearance.
    pub fn composition(&self) -> Vec<(Material, f64)> {
        let mut out: Vec<(Material, f64)> = Vec::new();
        for (k, l) in self.layers.iter().enumerate() {
            let m = l.material.density() * (sphere_volume(l.outer_radius) - sphere_volume(self.inner_radius(k)));
            if let Some(e) = out.iter_mut().find(|(mat, _)| *mat == l.material) {
                e.1 += m;
            } else {
                out.push((l.material.clone(), m));
            }
        }
        out
    }

    /// Mass-weighted mean coefficient of restitution.
    pub fn restitution(&self) -> f64 {
        let comp = self.composition();
        let total: f64 = comp.iter().map(|(_, m)| m).sum();
        comp.iter().map(|(mat, m)| mat.restitution() * m).sum::<f64>() / total
    }

    /// Homologous rescaling of all radii by `factor` (densities unchanged, so
    /// mass scales as `factor³`).
    pub fn scaled(&self, factor: f64) -> Self {
        SphericalStructure {
            layers: self
                .layers
                .iter()
                .map(|l| Layer { outer_radius: l.outer_radius * factor, material: l.material.clone() })
                .collect(),
        }
    }

    pub fn scaled_to_mass(&self, mass: f64) -> Result<Self, StructureError> {
        if !(mass > 0.0) || !mass.is_finite() {
            return Err(StructureError::BadMass(mass));
        }
        Ok(self.scaled((mass / self.mass()).cbrt()))
    }

    pub fn scaled_to_radius(&self, radius: f64) -> Result<Self, StructureError> {
        if !(radius > 0.0) || !radius.is_finite() {
            return Err(StructureError::BadRadius { index: self.layers.len() - 1, radius });
        }
        Ok(self.scaled(radius / self.radius()))
    }

    /// Replace the whole structure by a uniform sphere of `material` with the
    /// same total mass.
    pub fn with_material(&self, material: Material) -> Result<Self, StructureError> {
        Self::uniform_with_mass(material, self.mass())
    }

    /// Add a shell of the given thickness on the outside.
    pub fn push_shell(&mut self, thickness: f64, material: Material) -> Result<(), StructureError> {
        let r = self.radius() + thickness;
        if !(thickness > 0.0) || !thickness.is_finite() {
            return Err(StructureError::BadRadius { index: self.layers.len(), radius: r });
        }
        self.layers.push(Layer { outer_radius: r, material });
        Ok(())
    }

    /// Remove the outermost shell (keeping at least one layer).
    pub fn pop_shell(&mut self) -> Option<Layer> {
        if self.layers.len() > 1 {
            self.layers.pop()
        } else {
            None
        }
    }

    /// Fully differentiated merger: every material keeps its volume, and the
    /// shells are re-stacked densest-first so the heaviest material forms the
    /// new core. Mass and volume are both conserved exactly.
    pub fn merged_with(&self, other: &SphericalStructure) -> SphericalStructure {
        let mut buckets: Vec<(Material, f64)> = Vec::new();
        for s in [self, other] {
            for (k, l) in s.layers.iter().enumerate() {
                let v = sphere_volume(l.outer_radius) - sphere_volume(s.inner_radius(k));
                if let Some(e) = buckets.iter_mut().find(|(m, _)| *m == l.material) {
                    e.1 += v;
                } else {
                    buckets.push((l.material.clone(), v));
                }
            }
        }
        buckets.sort_by(|a, b| b.0.density().partial_cmp(&a.0.density()).unwrap_or(std::cmp::Ordering::Equal));
        let mut cumulative = 0.0;
        let layers = buckets
            .into_iter()
            .map(|(material, v)| {
                cumulative += v;
                Layer { outer_radius: radius_of_volume(cumulative), material }
            })
            .collect();
        SphericalStructure { layers }
    }
}

/// Two-layer (core + mantle) model reproducing a prescribed bulk mass and
/// radius. Solves `M = 4π/3 [ρ_c r_c³ + ρ_m (R³ − r_c³)]` for the core radius.
///
/// Falls back to a uniform sphere of a custom material with the exact bulk
/// density when no core fraction in `[0, 1]` fits (for instance Saturn is less
/// dense than any built-in material, and Uranus is less dense than `Gas`).
pub fn fit_two_layer(core: Material, mantle: Material, mass: f64, radius: f64) -> Result<SphericalStructure, StructureError> {
    if !(mass > 0.0) || !mass.is_finite() {
        return Err(StructureError::BadMass(mass));
    }
    if !(radius > 0.0) || !radius.is_finite() {
        return Err(StructureError::BadRadius { index: 0, radius });
    }
    let rho_c = core.density();
    let rho_m = mantle.density();
    let bulk = mass / sphere_volume(radius);
    let r3 = radius * radius * radius;
    if (rho_c - rho_m).abs() > 0.0 {
        let rc3 = (mass / (4.0 / 3.0 * PI) - rho_m * r3) / (rho_c - rho_m);
        if rc3 > 0.0 && rc3 < r3 {
            let rc = rc3.cbrt();
            return SphericalStructure::from_layers(vec![
                Layer { outer_radius: rc, material: core },
                Layer { outer_radius: radius, material: mantle },
            ]);
        }
    }
    let restitution = 0.5 * (core.restitution() + mantle.restitution());
    let mut custom = Material::custom(format!("{}/{}", core.name(), mantle.name()), bulk, restitution);
    if let Material::Custom { colour, .. } = &mut custom {
        *colour = mantle.colour();
    }
    SphericalStructure::uniform_with_radius(custom, radius)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{M_EARTH, R_EARTH};
    use approx::assert_relative_eq;

    #[test]
    fn uniform_sphere_mass_and_inertia() {
        let s = SphericalStructure::uniform_with_radius(Material::Rock, 1000.0).unwrap();
        assert_relative_eq!(s.mass(), 3000.0 * sphere_volume(1000.0), max_relative = 1e-12);
        assert_relative_eq!(s.inertia_factor(), 0.4, max_relative = 1e-12);
        assert_relative_eq!(s.mean_density(), 3000.0, max_relative = 1e-12);
    }

    #[test]
    fn uniform_from_mass_round_trips() {
        let s = SphericalStructure::uniform_with_mass(Material::Iron, 1e20).unwrap();
        assert_relative_eq!(s.mass(), 1e20, max_relative = 1e-12);
    }

    #[test]
    fn two_layer_fit_matches_bulk_properties() {
        let s = fit_two_layer(Material::Iron, Material::Rock, M_EARTH, R_EARTH).unwrap();
        assert_eq!(s.layers().len(), 2);
        assert_relative_eq!(s.mass(), M_EARTH, max_relative = 1e-10);
        assert_relative_eq!(s.radius(), R_EARTH, max_relative = 1e-12);
        // Centrally condensed: smaller inertia factor than a uniform sphere.
        assert!(s.inertia_factor() < 0.4);
        assert!(s.inertia_factor() > 0.3);
    }

    #[test]
    fn two_layer_fit_falls_back_when_too_light() {
        // Saturn: 687 kg/m^3 is lighter than any built-in material.
        let s = fit_two_layer(Material::Rock, Material::Gas, 5.6834e26, 5.8232e7).unwrap();
        assert_eq!(s.layers().len(), 1);
        assert_relative_eq!(s.mass(), 5.6834e26, max_relative = 1e-10);
    }

    #[test]
    fn enclosed_mass_is_monotone_and_saturates() {
        let s = fit_two_layer(Material::Iron, Material::Rock, M_EARTH, R_EARTH).unwrap();
        let mut prev = 0.0;
        for i in 1..=20 {
            let r = R_EARTH * i as f64 / 10.0;
            let m = s.enclosed_mass(r);
            assert!(m >= prev);
            prev = m;
        }
        assert_relative_eq!(prev, M_EARTH, max_relative = 1e-10);
    }

    #[test]
    fn merger_conserves_mass_and_volume_and_differentiates() {
        let a = SphericalStructure::uniform_with_radius(Material::Ice, 100.0).unwrap();
        let b = SphericalStructure::uniform_with_radius(Material::Iron, 50.0).unwrap();
        let m = a.merged_with(&b);
        assert_relative_eq!(m.mass(), a.mass() + b.mass(), max_relative = 1e-12);
        assert_relative_eq!(m.volume(), a.volume() + b.volume(), max_relative = 1e-12);
        assert_eq!(m.layers()[0].material, Material::Iron, "iron sinks to the core");
        assert_eq!(m.layers()[1].material, Material::Ice);
    }

    #[test]
    fn homologous_scaling() {
        let s = fit_two_layer(Material::Iron, Material::Rock, M_EARTH, R_EARTH).unwrap();
        let t = s.scaled_to_mass(8.0 * M_EARTH).unwrap();
        assert_relative_eq!(t.radius(), 2.0 * R_EARTH, max_relative = 1e-12);
        assert_relative_eq!(t.mass(), 8.0 * M_EARTH, max_relative = 1e-10);
    }

    #[test]
    fn invalid_layers_rejected() {
        assert!(SphericalStructure::from_layers(vec![]).is_err());
        assert!(SphericalStructure::from_layers(vec![
            Layer { outer_radius: 10.0, material: Material::Rock },
            Layer { outer_radius: 10.0, material: Material::Ice },
        ])
        .is_err());
        assert!(SphericalStructure::uniform_with_mass(Material::Rock, -1.0).is_err());
    }
}
