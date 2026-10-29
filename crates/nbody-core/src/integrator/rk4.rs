//! Classical fourth-order Runge–Kutta, included as the non-symplectic
//! control: for the same step size it is as accurate as Yoshida over one
//! orbit, but its energy error grows linearly in time, so orbits slowly
//! decay. Four force evaluations per step.

use super::{AccelFn, Integrator, IntegratorKind, State};
use crate::vector::Vec3;

#[derive(Debug, Default)]
pub struct Rk4 {
    a1: Vec<Vec3>,
    a2: Vec<Vec3>,
    a3: Vec<Vec3>,
    a4: Vec<Vec3>,
    v2: Vec<Vec3>,
    v3: Vec<Vec3>,
    v4: Vec<Vec3>,
    x: Vec<Vec3>,
}

impl Integrator for Rk4 {
    fn kind(&self) -> IntegratorKind {
        IntegratorKind::Rk4
    }

    fn step(&mut self, s: &mut State, dt: f64, accel: &mut AccelFn) {
        let n = s.len();
        s.ensure_acc(accel);
        self.a1.clear();
        self.a1.extend_from_slice(&s.acc);
        let h2 = 0.5 * dt;

        // Stage 2
        self.x.clear();
        self.v2.clear();
        for i in 0..n {
            self.x.push(s.pos[i] + s.vel[i] * h2);
            self.v2.push(s.vel[i] + self.a1[i] * h2);
        }
        accel(&self.x, &mut self.a2);

        // Stage 3
        self.x.clear();
        self.v3.clear();
        for i in 0..n {
            self.x.push(s.pos[i] + self.v2[i] * h2);
            self.v3.push(s.vel[i] + self.a2[i] * h2);
        }
        accel(&self.x, &mut self.a3);

        // Stage 4
        self.x.clear();
        self.v4.clear();
        for i in 0..n {
            self.x.push(s.pos[i] + self.v3[i] * dt);
            self.v4.push(s.vel[i] + self.a3[i] * dt);
        }
        accel(&self.x, &mut self.a4);

        let w = dt / 6.0;
        for i in 0..n {
            s.pos[i] += (s.vel[i] + self.v2[i] * 2.0 + self.v3[i] * 2.0 + self.v4[i]) * w;
            s.vel[i] += (self.a1[i] + self.a2[i] * 2.0 + self.a3[i] * 2.0 + self.a4[i]) * w;
        }
        s.acc_valid = false;
    }
}
