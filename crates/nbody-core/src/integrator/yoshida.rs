//! Yoshida's fourth-order symplectic integrator (Yoshida, *Phys. Lett. A*
//! 150, 262, 1990).
//!
//! The "triple jump": a second-order symmetric scheme `S₂(h)` composed as
//! `S₂(w₁h) S₂(w₀h) S₂(w₁h)` with
//!
//! ```text
//! w₁ = 1 / (2 − 2^{1/3}),   w₀ = −2^{1/3} / (2 − 2^{1/3}) = 1 − 2w₁
//! ```
//!
//! cancels the third-order error term and yields fourth order. The middle
//! sub-step runs *backwards in time* (w₀ < 0), which is characteristic of
//! all fourth-order splittings with positive-only stages being impossible.
//! Three force evaluations per step, the acceleration cache being reused
//! between consecutive kicks.

use super::{AccelFn, Integrator, IntegratorKind, State};

const CBRT2: f64 = 1.259_921_049_894_873_2;
pub const W1: f64 = 1.0 / (2.0 - CBRT2);
pub const W0: f64 = 1.0 - 2.0 * W1;

#[derive(Debug, Default, Clone, Copy)]
pub struct Yoshida4;

impl Integrator for Yoshida4 {
    fn kind(&self) -> IntegratorKind {
        IntegratorKind::Yoshida4
    }

    fn step(&mut self, s: &mut State, dt: f64, accel: &mut AccelFn) {
        // Three KDK leapfrogs; adjacent half-kicks are merged into one kick.
        s.ensure_acc(accel);
        s.kick(0.5 * W1 * dt);
        s.drift(W1 * dt);
        s.evaluate(accel);
        s.kick(0.5 * (W1 + W0) * dt);
        s.drift(W0 * dt);
        s.evaluate(accel);
        s.kick(0.5 * (W0 + W1) * dt);
        s.drift(W1 * dt);
        s.evaluate(accel);
        s.kick(0.5 * W1 * dt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weights_sum_to_one() {
        assert!((2.0 * W1 + W0 - 1.0).abs() < 1e-15);
        assert!((W1 - 1.351_207_191_959_657_6).abs() < 1e-12);
        assert!((W0 + 1.702_414_383_919_315_2).abs() < 1e-12);
    }
}
