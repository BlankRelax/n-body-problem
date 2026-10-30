//! The simulation driver: owns a system, an integrator, a force method, a
//! time-step policy and the collision handling, and keeps the bookkeeping
//! (energy error, collision log, projectile trackers) up to date.

use crate::collision::{self, CollisionEvent, CollisionMode};
use crate::diagnostics::{ConservationTracker, Diagnostics};
use crate::force::{self, ForceMethod};
use crate::integrator::{Integrator, IntegratorKind, State};
use crate::projectile::{Outcome, ProjectileError, ProjectileSpec, ProjectileTracker};
use crate::system::NBodySystem;
use crate::vector::Vec3;
use crate::body::BodyId;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TimestepPolicy {
    Fixed {
        dt: f64,
    },
    /// `dt = η · t_dyn` clamped to `[dt_min, dt_max]`, where `t_dyn` is the
    /// shortest pairwise two-body dynamical time. Breaks exact symplecticity
    /// (the step depends on the state) but keeps close encounters resolved.
    Adaptive {
        eta: f64,
        dt_min: f64,
        dt_max: f64,
    },
}

impl TimestepPolicy {
    pub fn fixed(dt: f64) -> Self {
        TimestepPolicy::Fixed { dt }
    }

    pub fn nominal_dt(&self) -> f64 {
        match *self {
            TimestepPolicy::Fixed { dt } => dt,
            TimestepPolicy::Adaptive { dt_max, .. } => dt_max,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SimulationConfig {
    pub integrator: IntegratorKind,
    pub timestep: TimestepPolicy,
    pub force: ForceMethod,
    pub collisions: CollisionMode,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        SimulationConfig {
            integrator: IntegratorKind::Yoshida4,
            timestep: TimestepPolicy::Fixed { dt: 1.0 },
            force: ForceMethod::Direct,
            collisions: CollisionMode::Physical,
        }
    }
}

pub struct Simulation {
    pub system: NBodySystem,
    config: SimulationConfig,
    integrator: Box<dyn Integrator>,
    state: State,
    mass: Vec<f64>,
    synced_generation: Option<u64>,
    pub steps: u64,
    pub force_evaluations: u64,
    pub last_dt: f64,
    pub collisions: Vec<CollisionEvent>,
    pub trackers: Vec<ProjectileTracker>,
    pub conservation: ConservationTracker,
}

impl std::fmt::Debug for Simulation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Simulation")
            .field("bodies", &self.system.len())
            .field("time", &self.system.time)
            .field("config", &self.config)
            .field("steps", &self.steps)
            .field("force_evaluations", &self.force_evaluations)
            .field("collisions", &self.collisions.len())
            .field("trackers", &self.trackers.len())
            .finish()
    }
}

impl Simulation {
    pub fn new(system: NBodySystem, config: SimulationConfig) -> Simulation {
        let conservation = ConservationTracker::new(&system);
        Simulation {
            integrator: config.integrator.build(),
            last_dt: config.timestep.nominal_dt(),
            system,
            config,
            state: State::default(),
            mass: Vec::new(),
            synced_generation: None,
            steps: 0,
            force_evaluations: 0,
            collisions: Vec::new(),
            trackers: Vec::new(),
            conservation,
        }
    }

    pub fn config(&self) -> &SimulationConfig {
        &self.config
    }

    pub fn set_config(&mut self, config: SimulationConfig) {
        if config.integrator != self.config.integrator {
            self.integrator = config.integrator.build();
            self.state.acc_valid = false;
        }
        self.config = config;
    }

    pub fn set_integrator(&mut self, kind: IntegratorKind) {
        let mut c = self.config.clone();
        c.integrator = kind;
        self.set_config(c);
    }

    pub fn set_timestep(&mut self, policy: TimestepPolicy) {
        self.config.timestep = policy;
    }

    pub fn set_collisions(&mut self, mode: CollisionMode) {
        self.config.collisions = mode;
    }

    pub fn set_force(&mut self, force: ForceMethod) {
        self.config.force = force;
        self.state.acc_valid = false;
    }

    /// Re-baseline the conservation tracker (after manual edits, say).
    pub fn reset_baseline(&mut self) {
        self.conservation = ConservationTracker::new(&self.system);
    }

    fn sync_state(&mut self) {
        if self.synced_generation != Some(self.system.generation()) || self.state.len() != self.system.len() {
            let (mut pos, mut vel) = (std::mem::take(&mut self.state.pos), std::mem::take(&mut self.state.vel));
            self.system.fill_state(&mut pos, &mut vel, &mut self.mass);
            self.state = State::new(pos, vel);
            self.synced_generation = Some(self.system.generation());
        }
    }

    /// Time step the policy selects for the current state.
    pub fn choose_dt(&mut self) -> f64 {
        self.sync_state();
        match self.config.timestep {
            TimestepPolicy::Fixed { dt } => dt,
            TimestepPolicy::Adaptive { eta, dt_min, dt_max } => {
                let t = force::min_dynamical_time(&self.state.pos, &self.mass, self.system.g, self.system.softening);
                if t.is_finite() {
                    (eta * t).clamp(dt_min, dt_max)
                } else {
                    dt_max
                }
            }
        }
    }

    /// Advance by one step; returns the step size used.
    pub fn step(&mut self) -> f64 {
        let dt = self.choose_dt();
        self.step_with(dt);
        dt
    }

    /// Advance by exactly `dt`.
    pub fn step_with(&mut self, dt: f64) {
        self.sync_state();
        if self.system.is_empty() {
            self.system.time += dt;
            return;
        }
        let g = self.system.g;
        let eps = self.system.softening;
        let force = self.config.force;
        let mass = &self.mass;
        let mut evals = 0u64;
        self.integrator.step(&mut self.state, dt, &mut |p: &[Vec3], out: &mut Vec<Vec3>| {
            force.accelerations(p, mass, g, eps, out);
            evals += 1;
        });
        self.force_evaluations += evals;
        self.system.write_state(&self.state.pos, &self.state.vel);
        self.system.time += dt;
        self.steps += 1;
        self.last_dt = dt;

        let mut new_events = Vec::new();
        if self.config.collisions != CollisionMode::None && !collision::detect_overlaps(self.system.bodies()).is_empty() {
            new_events = collision::resolve(&mut self.system, self.config.collisions);
            for e in &new_events {
                self.conservation.account_collision(e.energy_change);
            }
        }
        for t in &mut self.trackers {
            t.update(&self.system, &new_events);
        }
        self.collisions.extend(new_events);
    }

    /// Advance until `duration` has elapsed; returns the number of steps.
    pub fn advance(&mut self, duration: f64) -> u64 {
        let t_end = self.system.time + duration;
        self.run_until(t_end)
    }

    /// Advance until the clock reaches `t_end` (the final step is shortened
    /// so the end time is hit exactly).
    pub fn run_until(&mut self, t_end: f64) -> u64 {
        let mut n = 0;
        while self.system.time < t_end {
            let dt = self.choose_dt().min(t_end - self.system.time);
            if dt <= 0.0 {
                break;
            }
            self.step_with(dt);
            n += 1;
        }
        n
    }

    pub fn diagnostics(&self) -> Diagnostics {
        Diagnostics::compute(&self.system)
    }

    /// Relative energy error since the baseline, net of collision losses.
    pub fn energy_error(&mut self) -> f64 {
        let d = self.diagnostics();
        self.conservation.record(&d)
    }

    // ---- projectiles --------------------------------------------------------------

    /// Add a projectile on its incoming hyperbola and start tracking it.
    pub fn inject_projectile(&mut self, spec: &ProjectileSpec) -> Result<BodyId, ProjectileError> {
        let (id, tracker) = crate::projectile::inject(&mut self.system, spec)?;
        self.trackers.push(tracker);
        self.reset_baseline();
        Ok(id)
    }

    /// Add a projectile at an explicit position and velocity (a launch from
    /// inside the system) and start tracking it.
    pub fn inject_projectile_at(
        &mut self,
        name: impl Into<String>,
        structure: crate::structure::SphericalStructure,
        position: Vec3,
        velocity: Vec3,
    ) -> Result<BodyId, ProjectileError> {
        let (id, tracker) = crate::projectile::inject_at(&mut self.system, name, structure, position, velocity)?;
        self.trackers.push(tracker);
        self.reset_baseline();
        Ok(id)
    }

    pub fn tracker(&self, id: BodyId) -> Option<&ProjectileTracker> {
        self.trackers.iter().find(|t| t.id == id)
    }

    /// Latest projectile's outcome, if any projectile has been injected.
    pub fn projectile_outcome(&self) -> Option<&Outcome> {
        self.trackers.last().map(|t| &t.outcome)
    }

    /// Run until the latest projectile reaches a terminal outcome or
    /// `max_time` of simulation time has elapsed. Returns the outcome.
    pub fn run_projectile(&mut self, max_time: f64) -> Outcome {
        let t_end = self.system.time + max_time;
        while self.system.time < t_end {
            if let Some(t) = self.trackers.last() {
                if t.outcome.is_terminal() {
                    break;
                }
            } else {
                break;
            }
            let dt = self.choose_dt().min(t_end - self.system.time);
            if dt <= 0.0 {
                break;
            }
            self.step_with(dt);
        }
        self.trackers.last().map(|t| t.outcome.clone()).unwrap_or(Outcome::Undetermined)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::body::Body;
    use crate::material::Material;
    use crate::orbit;
    use crate::units::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    fn sun_earth() -> NBodySystem {
        let mut s = NBodySystem::si();
        let sun = s.add(Body::uniform("Sun", Material::Plasma, M_SUN, Vec3::ZERO, Vec3::ZERO).unwrap().with_kind(crate::body::BodyKind::Star));
        let earth = s.add(Body::uniform("Earth", Material::Rock, M_EARTH, Vec3::new(AU, 0.0, 0.0), Vec3::ZERO).unwrap());
        s.set_circular_orbit(earth, sun, Vec3::Z).unwrap();
        s.move_to_barycentre();
        s
    }

    #[test]
    fn earth_returns_after_one_year_with_tiny_energy_error() {
        let s = sun_earth();
        let x0 = s.find("Earth").unwrap().position;
        let mu = G_SI * (M_SUN + M_EARTH);
        let period = orbit::kepler_period(mu, AU);
        let cfg = SimulationConfig {
            integrator: IntegratorKind::Yoshida4,
            timestep: TimestepPolicy::fixed(period / 2000.0),
            force: ForceMethod::Direct,
            collisions: CollisionMode::Physical,
        };
        let mut sim = Simulation::new(s, cfg);
        sim.advance(period);
        let x1 = sim.system.find("Earth").unwrap().position;
        assert!((x1 - x0).norm() / AU < 1e-6, "position error {}", (x1 - x0).norm() / AU);
        assert!(sim.energy_error().abs() < 1e-10);
        assert!((sim.system.time - period).abs() < 1e-6);
    }

    #[test]
    fn adaptive_steps_shrink_near_periapsis() {
        let mut s = NBodySystem::new(1.0);
        s.add(Body::uniform("A", Material::Rock, 1.0, Vec3::ZERO, Vec3::ZERO).unwrap());
        // e = 0.9 orbit, start at apoapsis.
        let (r, v) = orbit::state_from_elements(
            1.0,
            &orbit::Elements { a: 1.0, e: 0.9, inclination: 0.0, raan: 0.0, arg_periapsis: 0.0, true_anomaly: PI },
        )
        .unwrap();
        s.add(Body::uniform("B", Material::Rock, 1e-9, r, v).unwrap());
        s.softening = 0.0;
        let cfg = SimulationConfig {
            integrator: IntegratorKind::Leapfrog,
            timestep: TimestepPolicy::Adaptive { eta: 0.02, dt_min: 1e-6, dt_max: 1.0 },
            force: ForceMethod::Direct,
            collisions: CollisionMode::None,
        };
        let mut sim = Simulation::new(s, cfg);
        let dt_apo = sim.choose_dt();
        sim.advance(PI); // half a period: now at periapsis
        let dt_peri = sim.choose_dt();
        assert!(dt_peri < dt_apo / 20.0, "apo {dt_apo} peri {dt_peri}");
        assert!(sim.energy_error().abs() < 1e-3);
    }

    #[test]
    fn collision_energy_is_booked_separately() {
        let mut s = NBodySystem::new(1.0);
        s.add(Body::uniform("A", Material::Rock, 1.0, Vec3::new(-1.0, 0.0, 0.0), Vec3::new(0.5, 0.0, 0.0)).unwrap());
        s.add(Body::uniform("B", Material::Rock, 1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(-0.5, 0.0, 0.0)).unwrap());
        let cfg = SimulationConfig {
            integrator: IntegratorKind::Pefrl,
            timestep: TimestepPolicy::Adaptive { eta: 0.005, dt_min: 1e-9, dt_max: 1e-2 },
            force: ForceMethod::Direct,
            collisions: CollisionMode::Merge,
        };
        let mut sim = Simulation::new(s, cfg);
        sim.advance(3.0);
        assert_eq!(sim.system.len(), 1);
        assert_eq!(sim.collisions.len(), 1);
        assert_relative_eq!(sim.conservation.collision_energy_change, sim.collisions[0].energy_change, max_relative = 1e-12);
        assert!(sim.conservation.collision_energy_change != 0.0);
        // Net of the booked collision change, the integration error is tiny.
        assert!(sim.energy_error().abs() < 1e-6, "numerical error {}", sim.energy_error());
    }
}
