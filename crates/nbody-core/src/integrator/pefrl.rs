//! Position-Extended Forest–Ruth-Like integrator (Omelyan, Mryglod & Folk,
//! *Comput. Phys. Commun.* 146, 188, 2002).
//!
//! A fourth-order symplectic splitting with four force evaluations per
//! step whose free coefficients were optimised to minimise the leading
//! error term; in practice it is ~30× more accurate than Yoshida's scheme
//! at equal step size, and still ~10× more accurate at equal cost.
//!
//! ```text
//! x += ξ h v
//! v += (1 − 2λ) h/2 a(x)
//! x += χ h v
//! v += λ h a(x)
//! x += (1 − 2(χ + ξ)) h v
//! v += λ h a(x)
//! x += χ h v
//! v += (1 − 2λ) h/2 a(x)
//! x += ξ h v
//! ```

use super::{AccelFn, Integrator, IntegratorKind, State};

pub const XI: f64 = 0.178_617_895_844_809_1;
pub const LAMBDA: f64 = -0.212_341_831_062_605_4;
pub const CHI: f64 = -0.066_264_582_669_818_49;

#[derive(Debug, Default, Clone, Copy)]
pub struct Pefrl;

impl Integrator for Pefrl {
    fn kind(&self) -> IntegratorKind {
        IntegratorKind::Pefrl
    }

    fn step(&mut self, s: &mut State, dt: f64, accel: &mut AccelFn) {
        s.drift(XI * dt);
        s.evaluate(accel);
        s.kick((1.0 - 2.0 * LAMBDA) * 0.5 * dt);
        s.drift(CHI * dt);
        s.evaluate(accel);
        s.kick(LAMBDA * dt);
        s.drift((1.0 - 2.0 * (CHI + XI)) * dt);
        s.evaluate(accel);
        s.kick(LAMBDA * dt);
        s.drift(CHI * dt);
        s.evaluate(accel);
        s.kick((1.0 - 2.0 * LAMBDA) * 0.5 * dt);
        s.drift(XI * dt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drift_coefficients_sum_to_one() {
        let total = 2.0 * XI + 2.0 * CHI + (1.0 - 2.0 * (CHI + XI));
        assert!((total - 1.0).abs() < 1e-15);
        let kicks = 2.0 * (1.0 - 2.0 * LAMBDA) * 0.5 + 2.0 * LAMBDA;
        assert!((kicks - 1.0).abs() < 1e-15);
    }
}
