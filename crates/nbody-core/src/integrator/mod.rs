//! Time integrators for the Hamiltonian system `H = Σ p²/2m + U(q)`.
//!
//! Gravity is a separable Hamiltonian, so it can be advanced with
//! *splitting* methods that alternate exact "kick" (momentum) and "drift"
//! (position) maps. Such compositions are **symplectic**: they conserve a
//! nearby shadow Hamiltonian exactly, so the energy error stays bounded
//! for exponentially long times instead of growing secularly. That is why
//! the symplectic schemes here are the right default for orbital dynamics,
//! while a classical Runge–Kutta method (also provided, for comparison)
//! slowly spirals orbits inwards.
//!
//! | Integrator | Order | Symplectic | Force evaluations / step |
//! |------------|-------|------------|--------------------------|
//! | Leapfrog (kick–drift–kick) | 2 | yes | 1 |
//! | Yoshida 4th order          | 4 | yes | 3 |
//! | PEFRL                      | 4 | yes | 4 (best error constant) |
//! | RK4                        | 4 | no  | 4 |
//!
//! Every scheme keeps the acceleration of the *current* positions cached in
//! `State::acc` when it is cheap to do so, so that a kick-first scheme does
//! not need an extra force evaluation at the start of each step.

pub mod leapfrog;
pub mod pefrl;
pub mod rk4;
pub mod yoshida;

use crate::vector::Vec3;
use serde::{Deserialize, Serialize};

/// Callback computing accelerations for a set of positions.
pub type AccelFn<'a> = dyn FnMut(&[Vec3], &mut Vec<Vec3>) + 'a;

/// Structure-of-arrays phase-space state used inside integrators.
#[derive(Debug, Clone, Default)]
pub struct State {
    pub pos: Vec<Vec3>,
    pub vel: Vec<Vec3>,
    pub acc: Vec<Vec3>,
    /// Whether `acc` corresponds to `pos`.
    pub acc_valid: bool,
}

impl State {
    pub fn new(pos: Vec<Vec3>, vel: Vec<Vec3>) -> State {
        let n = pos.len();
        State { pos, vel, acc: vec![Vec3::ZERO; n], acc_valid: false }
    }

    pub fn len(&self) -> usize {
        self.pos.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }

    /// Make sure `acc` matches `pos`, evaluating the force if needed.
    pub fn ensure_acc(&mut self, accel: &mut AccelFn) {
        if !self.acc_valid || self.acc.len() != self.pos.len() {
            accel(&self.pos, &mut self.acc);
            self.acc_valid = true;
        }
    }

    #[inline]
    pub fn kick(&mut self, dt: f64) {
        for (v, a) in self.vel.iter_mut().zip(&self.acc) {
            *v += *a * dt;
        }
    }

    #[inline]
    pub fn drift(&mut self, dt: f64) {
        for (p, v) in self.pos.iter_mut().zip(&self.vel) {
            *p += *v * dt;
        }
        self.acc_valid = false;
    }

    /// Evaluate the force at the current positions.
    #[inline]
    pub fn evaluate(&mut self, accel: &mut AccelFn) {
        accel(&self.pos, &mut self.acc);
        self.acc_valid = true;
    }
}

pub trait Integrator: Send {
    fn kind(&self) -> IntegratorKind;
    /// Advance the state by `dt`.
    fn step(&mut self, s: &mut State, dt: f64, accel: &mut AccelFn);
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum IntegratorKind {
    Leapfrog,
    #[default]
    Yoshida4,
    Pefrl,
    Rk4,
}

impl IntegratorKind {
    pub const ALL: [IntegratorKind; 4] = [IntegratorKind::Leapfrog, IntegratorKind::Yoshida4, IntegratorKind::Pefrl, IntegratorKind::Rk4];

    pub fn name(self) -> &'static str {
        match self {
            IntegratorKind::Leapfrog => "leapfrog",
            IntegratorKind::Yoshida4 => "yoshida4",
            IntegratorKind::Pefrl => "pefrl",
            IntegratorKind::Rk4 => "rk4",
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            IntegratorKind::Leapfrog => "kick-drift-kick leapfrog: 2nd order, symplectic, 1 force evaluation",
            IntegratorKind::Yoshida4 => "Yoshida (1990) triple-jump composition: 4th order, symplectic, 3 force evaluations",
            IntegratorKind::Pefrl => "Omelyan-Mryglod-Folk (2002) PEFRL: 4th order, symplectic, 4 force evaluations, smallest error constant",
            IntegratorKind::Rk4 => "classical Runge-Kutta: 4th order, NOT symplectic (energy drifts), 4 force evaluations",
        }
    }

    pub fn order(self) -> u32 {
        match self {
            IntegratorKind::Leapfrog => 2,
            _ => 4,
        }
    }

    pub fn symplectic(self) -> bool {
        !matches!(self, IntegratorKind::Rk4)
    }

    pub fn force_evaluations_per_step(self) -> u32 {
        match self {
            IntegratorKind::Leapfrog => 1,
            IntegratorKind::Yoshida4 => 3,
            IntegratorKind::Pefrl => 4,
            IntegratorKind::Rk4 => 4,
        }
    }

    pub fn parse(s: &str) -> Option<IntegratorKind> {
        Some(match s.trim().to_ascii_lowercase().as_str() {
            "leapfrog" | "kdk" | "verlet" | "lf" => IntegratorKind::Leapfrog,
            "yoshida" | "yoshida4" | "y4" => IntegratorKind::Yoshida4,
            "pefrl" | "forest-ruth" | "fr" => IntegratorKind::Pefrl,
            "rk4" | "runge-kutta" => IntegratorKind::Rk4,
            _ => return None,
        })
    }

    pub fn build(self) -> Box<dyn Integrator> {
        match self {
            IntegratorKind::Leapfrog => Box::new(leapfrog::Leapfrog),
            IntegratorKind::Yoshida4 => Box::new(yoshida::Yoshida4),
            IntegratorKind::Pefrl => Box::new(pefrl::Pefrl),
            IntegratorKind::Rk4 => Box::new(rk4::Rk4::default()),
        }
    }
}

impl std::fmt::Display for IntegratorKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Kepler problem with μ = 1: a unit circular orbit has period 2π.
    pub(crate) fn kepler_accel(pos: &[Vec3], out: &mut Vec<Vec3>) {
        out.clear();
        for p in pos {
            let r2 = p.norm_sq();
            out.push(-*p / (r2 * r2.sqrt()));
        }
    }

    fn energy(s: &State) -> f64 {
        0.5 * s.vel[0].norm_sq() - 1.0 / s.pos[0].norm()
    }

    fn run(kind: IntegratorKind, dt: f64, orbits: f64, e: f64) -> (f64, f64) {
        // Start at periapsis of an orbit with semi-major axis 1 and eccentricity e.
        let rp = 1.0 - e;
        let vp = ((1.0 + e) / rp).sqrt();
        let mut s = State::new(vec![Vec3::new(rp, 0.0, 0.0)], vec![Vec3::new(0.0, vp, 0.0)]);
        let mut integ = kind.build();
        let e0 = energy(&s);
        let steps = (orbits * 2.0 * PI / dt).round() as usize;
        let mut max_drift: f64 = 0.0;
        for _ in 0..steps {
            integ.step(&mut s, dt, &mut kepler_accel);
            max_drift = max_drift.max(((energy(&s) - e0) / e0).abs());
        }
        let final_err = (s.pos[0] - Vec3::new(rp, 0.0, 0.0)).norm();
        (max_drift, final_err)
    }

    #[test]
    fn all_integrators_close_a_circular_orbit() {
        for kind in IntegratorKind::ALL {
            let (drift, err) = run(kind, 2.0 * PI / 400.0, 10.0, 0.0);
            let tol = if kind.order() == 2 { 1e-3 } else { 1e-6 };
            assert!(drift < tol, "{kind}: energy drift {drift}");
            assert!(err < 10.0 * tol.sqrt(), "{kind}: position error {err}");
        }
    }

    #[test]
    fn symplectic_energy_error_is_bounded_but_rk4_drifts() {
        // 200 orbits of an e=0.5 ellipse with a coarse step.
        let dt = 2.0 * PI / 300.0;
        let (lf, _) = run(IntegratorKind::Leapfrog, dt, 200.0, 0.5);
        let (y4, _) = run(IntegratorKind::Yoshida4, dt, 200.0, 0.5);
        let (pe, _) = run(IntegratorKind::Pefrl, dt, 200.0, 0.5);
        let (rk, _) = run(IntegratorKind::Rk4, dt, 200.0, 0.5);
        assert!(lf < 2e-2, "leapfrog {lf}");
        assert!(y4 < 1e-4, "yoshida {y4}");
        assert!(pe < y4, "pefrl ({pe}) should beat yoshida ({y4})");
        assert!(rk > y4, "rk4 ({rk}) should drift more than yoshida ({y4})");
    }

    #[test]
    fn convergence_orders() {
        for kind in IntegratorKind::ALL {
            let dt1 = 2.0 * PI / 100.0;
            let (_, e1) = run(kind, dt1, 1.0, 0.3);
            let (_, e2) = run(kind, dt1 / 2.0, 1.0, 0.3);
            let observed = (e1 / e2).log2();
            assert!(observed > kind.order() as f64 - 0.5, "{kind}: observed order {observed}");
        }
    }
}
