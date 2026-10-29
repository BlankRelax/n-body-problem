//! Gravitational acceleration kernels.
//!
//! Two evaluators share one physical model, the Plummer-softened pair force
//!
//! ```text
//! a_i = −G Σ_j m_j (r_i − r_j) / (|r_i − r_j|² + ε²)^{3/2}
//! ```
//!
//! * [`direct`] sums every pair exactly: O(N²), the right choice below a few
//!   thousand bodies and the reference the tree is tested against.
//! * [`barnes_hut`] groups distant bodies into multipole (monopole) nodes of
//!   an octree: O(N log N) with a controllable error set by the opening angle.
//!
//! Softening removes the singularity of the point-mass potential; with ε = 0
//! the kernel is exact Newtonian gravity. Collisions between finite-size
//! bodies are handled separately (see [`crate::collision`]), so ε can stay
//! zero for planetary systems and is mainly useful for collisionless clusters.

pub mod barnes_hut;
pub mod direct;

use crate::vector::Vec3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ForceMethod {
    Direct,
    /// Octree with multipole-acceptance criterion `s/d < θ`.
    BarnesHut { theta: f64 },
}

impl Default for ForceMethod {
    fn default() -> Self {
        ForceMethod::Direct
    }
}

impl ForceMethod {
    pub fn name(&self) -> &'static str {
        match self {
            ForceMethod::Direct => "direct",
            ForceMethod::BarnesHut { .. } => "barnes_hut",
        }
    }

    /// Fill `out` with the acceleration of every body.
    pub fn accelerations(&self, pos: &[Vec3], mass: &[f64], g: f64, softening: f64, out: &mut Vec<Vec3>) {
        match *self {
            ForceMethod::Direct => direct::accelerations(pos, mass, g, softening, out),
            ForceMethod::BarnesHut { theta } => barnes_hut::accelerations(pos, mass, g, softening, theta, out),
        }
    }
}

/// Plummer-softened potential energy `−G Σ_{i<j} m_i m_j / sqrt(r² + ε²)`,
/// the exact conjugate of the softened force so that `K + U` is conserved.
pub fn potential_energy(pos: &[Vec3], mass: &[f64], g: f64, softening: f64) -> f64 {
    let n = pos.len();
    let eps2 = softening * softening;
    let mut u = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let r2 = pos[i].distance_sq(pos[j]) + eps2;
            if r2 > 0.0 {
                u -= g * mass[i] * mass[j] / r2.sqrt();
            }
        }
    }
    u
}

/// Shortest two-body dynamical time over all pairs,
/// `min_{i<j} sqrt((r_ij² + ε²)^{3/2} / (G (m_i + m_j)))`, i.e. `1/ω` of the
/// tightest local Kepler problem. Used by the adaptive time-step policy.
pub fn min_dynamical_time(pos: &[Vec3], mass: &[f64], g: f64, softening: f64) -> f64 {
    let n = pos.len();
    let eps2 = softening * softening;
    let mut best = f64::INFINITY;
    for i in 0..n {
        for j in (i + 1)..n {
            let r2 = pos[i].distance_sq(pos[j]) + eps2;
            let msum = mass[i] + mass[j];
            if r2 > 0.0 && msum > 0.0 {
                let t2 = r2 * r2.sqrt() / (g * msum);
                if t2 < best {
                    best = t2;
                }
            }
        }
    }
    best.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The force must be minus the gradient of the potential; check by central
    /// finite differences on a random 5-body configuration with softening.
    #[test]
    fn force_is_gradient_of_potential() {
        let pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.2, -0.3),
            Vec3::new(-0.7, 1.1, 0.4),
            Vec3::new(0.3, -0.9, 1.2),
            Vec3::new(-1.3, -0.4, -0.8),
        ];
        let mass = vec![3.0, 1.0, 0.5, 2.0, 0.7];
        let (g, eps) = (1.3, 0.1);
        let mut acc = Vec::new();
        direct::accelerations(&pos, &mass, g, eps, &mut acc);
        let h = 1e-6;
        for i in 0..pos.len() {
            for k in 0..3 {
                let mut p = pos.clone();
                p[i][k] += h;
                let up = potential_energy(&p, &mass, g, eps);
                p[i][k] -= 2.0 * h;
                let um = potential_energy(&p, &mass, g, eps);
                let f_fd = -(up - um) / (2.0 * h);
                let f = mass[i] * acc[i][k];
                assert!((f - f_fd).abs() < 1e-6 * f.abs().max(1.0), "body {i} axis {k}: {f} vs {f_fd}");
            }
        }
    }

    #[test]
    fn dynamical_time_of_unit_binary() {
        let pos = vec![Vec3::ZERO, Vec3::X];
        let mass = vec![0.5, 0.5];
        assert!((min_dynamical_time(&pos, &mass, 1.0, 0.0) - 1.0).abs() < 1e-12);
    }
}
