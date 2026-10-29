//! Conserved quantities and integration-quality metrics.
//!
//! For an isolated Newtonian system the total energy, linear momentum and
//! angular momentum are exact invariants, so their drift measures the
//! integrator's error. Collisions change the energy *physically*
//! (inelastic mergers dissipate it); the tracker books those changes
//! separately so the reported error is purely numerical.

use crate::system::NBodySystem;
use crate::vector::Vec3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Diagnostics {
    pub time: f64,
    pub bodies: usize,
    pub total_mass: f64,
    pub kinetic: f64,
    pub potential: f64,
    pub total: f64,
    pub momentum: Vec3,
    pub angular_momentum: Vec3,
    pub center_of_mass: Vec3,
    pub com_velocity: Vec3,
    /// `2K / |U|`: 1 for a system in virial equilibrium, >1 if it is flying apart.
    pub virial_ratio: f64,
}

impl Diagnostics {
    pub fn compute(s: &NBodySystem) -> Diagnostics {
        let kinetic = s.kinetic_energy();
        let potential = s.potential_energy();
        Diagnostics {
            time: s.time,
            bodies: s.len(),
            total_mass: s.total_mass(),
            kinetic,
            potential,
            total: kinetic + potential,
            momentum: s.momentum(),
            angular_momentum: s.angular_momentum(),
            center_of_mass: s.center_of_mass(),
            com_velocity: s.com_velocity(),
            virial_ratio: if potential != 0.0 { 2.0 * kinetic / potential.abs() } else { f64::NAN },
        }
    }
}

/// Tracks drift of the invariants relative to a baseline, net of collisions.
#[derive(Debug, Clone)]
pub struct ConservationTracker {
    pub initial: Diagnostics,
    /// Sum of the energy changes caused by collision events since the baseline.
    pub collision_energy_change: f64,
    pub max_energy_error: f64,
    pub max_momentum_error: f64,
    pub max_angular_momentum_error: f64,
    energy_scale: f64,
    momentum_scale: f64,
    angular_scale: f64,
}

impl ConservationTracker {
    pub fn new(s: &NBodySystem) -> ConservationTracker {
        let initial = Diagnostics::compute(s);
        let energy_scale = initial.total.abs().max(initial.kinetic.abs()).max(initial.potential.abs()).max(f64::MIN_POSITIVE);
        // Characteristic momentum: the larger of the actual Σ m|v| and the
        // momentum a body would have with the system's energy scale, so that a
        // system starting at rest still gets a finite, meaningful scale.
        let extent = s.extent().max(f64::MIN_POSITIVE);
        let momentum_scale = s
            .bodies()
            .iter()
            .map(|b| b.mass() * b.velocity.norm())
            .sum::<f64>()
            .max((2.0 * initial.total_mass * energy_scale).sqrt())
            .max(f64::MIN_POSITIVE);
        let angular_scale = s
            .bodies()
            .iter()
            .map(|b| b.mass() * b.position.norm() * b.velocity.norm())
            .sum::<f64>()
            .max(momentum_scale * extent)
            .max(f64::MIN_POSITIVE);
        ConservationTracker {
            initial,
            collision_energy_change: 0.0,
            max_energy_error: 0.0,
            max_momentum_error: 0.0,
            max_angular_momentum_error: 0.0,
            energy_scale,
            momentum_scale,
            angular_scale,
        }
    }

    pub fn account_collision(&mut self, energy_change: f64) {
        self.collision_energy_change += energy_change;
    }

    /// Relative energy error `(E − E₀ − ΔE_collisions) / |E₀|`.
    pub fn energy_error(&self, now: &Diagnostics) -> f64 {
        (now.total - self.initial.total - self.collision_energy_change) / self.energy_scale
    }

    pub fn momentum_error(&self, now: &Diagnostics) -> f64 {
        (now.momentum - self.initial.momentum).norm() / self.momentum_scale
    }

    pub fn angular_momentum_error(&self, now: &Diagnostics) -> f64 {
        (now.angular_momentum - self.initial.angular_momentum).norm() / self.angular_scale
    }

    /// Update the running maxima; returns the current relative energy error.
    pub fn record(&mut self, now: &Diagnostics) -> f64 {
        let e = self.energy_error(now);
        self.max_energy_error = self.max_energy_error.max(e.abs());
        self.max_momentum_error = self.max_momentum_error.max(self.momentum_error(now));
        self.max_angular_momentum_error = self.max_angular_momentum_error.max(self.angular_momentum_error(now));
        e
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::body::Body;
    use crate::material::Material;

    #[test]
    fn virial_ratio_of_circular_binary_is_one() {
        let mut s = NBodySystem::new(1.0);
        let a = s.add(Body::uniform("A", Material::Rock, 1.0, Vec3::ZERO, Vec3::ZERO).unwrap());
        let b = s.add(Body::uniform("B", Material::Rock, 1e-6, Vec3::X, Vec3::ZERO).unwrap());
        s.set_circular_orbit(b, a, Vec3::Z).unwrap();
        s.move_to_barycentre();
        let d = Diagnostics::compute(&s);
        assert!((d.virial_ratio - 1.0).abs() < 1e-5, "virial {}", d.virial_ratio);
        assert!(d.total < 0.0);
    }
}
