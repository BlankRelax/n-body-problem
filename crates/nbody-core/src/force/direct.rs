//! Exact pairwise summation, O(N²).
//!
//! Each pair is visited once and Newton's third law is applied, halving the
//! work; the inner loop is branch-free apart from the coincident-point guard.

use crate::vector::Vec3;

/// Softened accelerations of every body by direct summation.
pub fn accelerations(pos: &[Vec3], mass: &[f64], g: f64, softening: f64, out: &mut Vec<Vec3>) {
    let n = pos.len();
    out.clear();
    out.resize(n, Vec3::ZERO);
    let eps2 = softening * softening;
    for i in 0..n {
        let pi = pos[i];
        let mut ai = Vec3::ZERO;
        for j in (i + 1)..n {
            let d = pos[j] - pi;
            let r2 = d.norm_sq() + eps2;
            if r2 <= 0.0 {
                continue;
            }
            let inv_r3 = 1.0 / (r2 * r2.sqrt());
            ai += d * (g * mass[j] * inv_r3);
            out[j] -= d * (g * mass[i] * inv_r3);
        }
        out[i] += ai;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_bodies_attract_with_inverse_square() {
        let pos = vec![Vec3::ZERO, Vec3::new(2.0, 0.0, 0.0)];
        let mass = vec![3.0, 5.0];
        let mut acc = Vec::new();
        accelerations(&pos, &mass, 1.0, 0.0, &mut acc);
        assert!((acc[0].x - 5.0 / 4.0).abs() < 1e-15);
        assert!((acc[1].x + 3.0 / 4.0).abs() < 1e-15);
        // Momentum conservation: Σ m a = 0.
        let p: Vec3 = acc.iter().zip(&mass).map(|(a, m)| *a * *m).sum();
        assert!(p.norm() < 1e-15);
    }

    #[test]
    fn coincident_bodies_do_not_produce_nan() {
        let pos = vec![Vec3::ZERO, Vec3::ZERO];
        let mass = vec![1.0, 1.0];
        let mut acc = Vec::new();
        accelerations(&pos, &mass, 1.0, 0.0, &mut acc);
        assert!(acc[0].is_finite() && acc[1].is_finite());
    }
}
