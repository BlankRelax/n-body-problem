//! Kick–drift–kick leapfrog (velocity Verlet).
//!
//! ```text
//! v_{n+1/2} = v_n + a(x_n) Δt/2
//! x_{n+1}   = x_n + v_{n+1/2} Δt
//! v_{n+1}   = v_{n+1/2} + a(x_{n+1}) Δt/2
//! ```
//!
//! Second order, symplectic and time-reversible. One force evaluation per
//! step because `a(x_{n+1})` is cached for the next step's first kick.

use super::{AccelFn, Integrator, IntegratorKind, State};

#[derive(Debug, Default, Clone, Copy)]
pub struct Leapfrog;

impl Integrator for Leapfrog {
    fn kind(&self) -> IntegratorKind {
        IntegratorKind::Leapfrog
    }

    fn step(&mut self, s: &mut State, dt: f64, accel: &mut AccelFn) {
        s.ensure_acc(accel);
        s.kick(0.5 * dt);
        s.drift(dt);
        s.evaluate(accel);
        s.kick(0.5 * dt);
    }
}
