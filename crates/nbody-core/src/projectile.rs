//! Projectile mode: inject an interloper on a hyperbolic approach and
//! classify what the system does to it.
//!
//! # Setting up the approach
//!
//! An object arriving from infinity is characterised by its excess speed
//! `v_∞` and impact parameter `b` (the offset of its undeflected straight
//! path from the system's centre of mass). Relative to the barycentre and
//! total mass of the system the incoming leg is a two-body hyperbola with
//!
//! ```text
//! a = μ / v_∞²,   e = sqrt(1 + b²/a²),   p = b²/a,   q = sqrt(a² + b²) − a
//! ```
//!
//! [`hyperbolic_approach_state`] evaluates that hyperbola at a chosen start
//! distance, so the projectile begins exactly on the trajectory it would
//! have followed from infinity instead of being launched with an ad-hoc
//! velocity. The start distance defaults to ten times the system's extent,
//! where the system's field is a monopole to better than a percent.
//!
//! # Classifying the outcome
//!
//! Energy conservation forbids capture by a *single* point mass: a body
//! arriving from infinity leaves again with the same `v_∞`. Capture needs a
//! third body to absorb energy (a close pass behind a planet is the classic
//! mechanism), and the same exchange can kick a planet out instead. The
//! [`ProjectileTracker`] therefore watches the projectile's specific energy
//! relative to the rest of the system,
//!
//! ```text
//! ε = ½ |v − v_cm|² − Σ_j G m_j / |r − r_j|
//! ```
//!
//! and reports
//!
//! * **Ejected** once ε > 0, the projectile is receding and has gone back
//!   past its start distance (with the outgoing `v_∞` and deflection angle);
//! * **Captured** once ε < 0 and the projectile has turned around at an
//!   apoapsis, i.e. it has started a second orbit;
//! * **Merged** when a collision absorbed it into another body;
//! * which original members, if any, were themselves unbound by the visit.

use crate::body::{Body, BodyId, BodyKind};
use crate::collision::{CollisionEvent, CollisionOutcome};
use crate::orbit::{self, OrbitalElements};
use crate::simulation::{Simulation, SimulationConfig};
use crate::structure::{SphericalStructure, StructureError};
use crate::system::NBodySystem;
use crate::vector::Vec3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ProjectileError {
    #[error("excess speed must be positive (got {0})")]
    NonPositiveSpeed(f64),
    #[error("impact parameter must be non-negative (got {0})")]
    NegativeImpactParameter(f64),
    #[error("start distance {start:.4e} m is inside the periapsis {periapsis:.4e} m of the approach hyperbola")]
    StartInsidePeriapsis { start: f64, periapsis: f64 },
    #[error("the target system has no bodies")]
    EmptySystem,
    #[error("approach direction must be a non-zero vector")]
    BadDirection,
    #[error(transparent)]
    Structure(#[from] StructureError),
}

#[derive(Debug, Clone, PartialEq)]
pub struct ProjectileSpec {
    pub name: String,
    pub structure: SphericalStructure,
    /// Hyperbolic excess speed, m/s.
    pub v_infinity: f64,
    /// Impact parameter, m.
    pub impact_parameter: f64,
    /// Direction of the velocity at infinity (need not be normalised).
    pub direction: Vec3,
    /// Direction of the impact-parameter offset; defaults to `z × direction`
    /// (in-plane for a system orbiting in the x-y plane).
    pub offset_direction: Option<Vec3>,
    /// Distance from the barycentre at which the projectile is created.
    pub start_distance: Option<f64>,
}

impl ProjectileSpec {
    pub fn new(name: impl Into<String>, structure: SphericalStructure, v_infinity: f64, impact_parameter: f64) -> Self {
        ProjectileSpec {
            name: name.into(),
            structure,
            v_infinity,
            impact_parameter,
            direction: Vec3::X,
            offset_direction: None,
            start_distance: None,
        }
    }

    pub fn with_direction(mut self, d: Vec3) -> Self {
        self.direction = d;
        self
    }

    pub fn with_offset_direction(mut self, d: Vec3) -> Self {
        self.offset_direction = Some(d);
        self
    }

    pub fn with_start_distance(mut self, d: f64) -> Self {
        self.start_distance = Some(d);
        self
    }
}

/// The incoming hyperbola in its perifocal frame plus the rotation into the
/// world frame. `None` for the degenerate radial case (`b = 0`).
struct Hyperbola {
    a: f64,
    e: f64,
    p: f64,
    q: f64,
    /// World-frame images of the perifocal x, y, z axes.
    axes: [Vec3; 3],
}

impl Hyperbola {
    fn rotate(&self, u: Vec3) -> Vec3 {
        self.axes[0] * u.x + self.axes[1] * u.y + self.axes[2] * u.z
    }
}

fn approach_directions(direction: Vec3, offset: Option<Vec3>) -> Result<(Vec3, Vec3), ProjectileError> {
    let dir = direction.try_normalized().ok_or(ProjectileError::BadDirection)?;
    let off = offset
        .unwrap_or_else(|| Vec3::Z.cross(dir))
        .reject_from_unit(dir)
        .try_normalized()
        .unwrap_or_else(|| dir.any_orthonormal());
    Ok((dir, off))
}

fn hyperbola(mu: f64, v_inf: f64, b: f64, dir: Vec3, off: Vec3) -> Option<Hyperbola> {
    let a = mu / (v_inf * v_inf);
    if b * v_inf * v_inf <= 1e-12 * mu.max(f64::MIN_POSITIVE) || b == 0.0 {
        return None;
    }
    let e = (1.0 + (b / a) * (b / a)).sqrt();
    let p = b * b / a;
    let q = p / (1.0 + e);
    // Incoming asymptote direction in the perifocal frame and the offset of
    // the asymptote from the focus (r_perp = v̂ × h / |v| with h along +z).
    let v_inf_pf = Vec3::new(1.0 / e, (e * e - 1.0).sqrt() / e, 0.0);
    let b_pf = Vec3::new(v_inf_pf.y, -v_inf_pf.x, 0.0);
    let n_pf = v_inf_pf.cross(b_pf);
    let n_world = dir.cross(off);
    let rotate = |u: Vec3| dir * u.dot(v_inf_pf) + off * u.dot(b_pf) + n_world * u.dot(n_pf);
    Some(Hyperbola { a, e, p, q, axes: [rotate(Vec3::X), rotate(Vec3::Y), rotate(Vec3::Z)] })
}

/// Position and velocity (relative to the attracting centre) of a body on
/// the incoming branch of the hyperbola defined by `v_inf` and `b`, at
/// distance `d`. `direction` is the velocity direction at infinity and
/// `offset` the direction in which the undeflected path is displaced.
pub fn hyperbolic_approach_state(
    mu: f64,
    v_inf: f64,
    b: f64,
    d: f64,
    direction: Vec3,
    offset: Option<Vec3>,
) -> Result<(Vec3, Vec3), ProjectileError> {
    Ok(approach_geometry(mu, v_inf, b, d, direction, offset)?.state())
}

/// Everything about an incoming hyperbola that is fixed analytically: the
/// start state, where and when periapsis is reached, and the deflection.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ApproachGeometry {
    pub position: Vec3,
    pub velocity: Vec3,
    /// Periapsis point relative to the focus.
    pub periapsis: Vec3,
    pub periapsis_distance: f64,
    /// Flight time from the start state to periapsis.
    pub time_to_periapsis: f64,
    /// Two-body deflection angle, rad.
    pub deflection: f64,
    pub eccentricity: f64,
}

impl ApproachGeometry {
    pub fn state(&self) -> (Vec3, Vec3) {
        (self.position, self.velocity)
    }
}

pub fn approach_geometry(
    mu: f64,
    v_inf: f64,
    b: f64,
    d: f64,
    direction: Vec3,
    offset: Option<Vec3>,
) -> Result<ApproachGeometry, ProjectileError> {
    if !(v_inf > 0.0) {
        return Err(ProjectileError::NonPositiveSpeed(v_inf));
    }
    if !(b >= 0.0) {
        return Err(ProjectileError::NegativeImpactParameter(b));
    }
    let (dir, off) = approach_directions(direction, offset)?;
    let Some(h) = hyperbola(mu, v_inf, b, dir, off) else {
        // Radial (degenerate) hyperbola: energy conservation fixes the speed;
        // the infall time is ∫ dr / sqrt(v_inf² + 2μ/r), done by Simpson's rule.
        let speed = (v_inf * v_inf + 2.0 * mu / d).sqrt();
        let n = 2000;
        let hstep = d / n as f64;
        let f = |r: f64| if r <= 0.0 { 0.0 } else { 1.0 / (v_inf * v_inf + 2.0 * mu / r).sqrt() };
        let mut t = f(0.0) + f(d);
        for k in 1..n {
            t += f(k as f64 * hstep) * if k % 2 == 1 { 4.0 } else { 2.0 };
        }
        t *= hstep / 3.0;
        return Ok(ApproachGeometry {
            position: -dir * d,
            velocity: dir * speed,
            periapsis: Vec3::ZERO,
            periapsis_distance: 0.0,
            time_to_periapsis: t,
            deflection: std::f64::consts::PI,
            eccentricity: 1.0,
        });
    };
    if d < h.q {
        return Err(ProjectileError::StartInsidePeriapsis { start: d, periapsis: h.q });
    }
    let cos_nu = ((h.p / d - 1.0) / h.e).clamp(-1.0, 1.0);
    let nu = -cos_nu.acos();
    let r_pf = Vec3::new(d * nu.cos(), d * nu.sin(), 0.0);
    let k = (mu / h.p).sqrt();
    let v_pf = Vec3::new(-k * nu.sin(), k * (h.e + nu.cos()), 0.0);
    // Hyperbolic Kepler equation gives the time from the start to periapsis.
    let big_h = 2.0 * (((h.e - 1.0) / (h.e + 1.0)).sqrt() * (0.5 * nu).tan()).atanh();
    let mean = h.e * big_h.sinh() - big_h;
    let time_to_periapsis = -mean * (h.a * h.a * h.a / mu).sqrt();
    Ok(ApproachGeometry {
        position: h.rotate(r_pf),
        velocity: h.rotate(v_pf),
        periapsis: h.rotate(Vec3::new(h.q, 0.0, 0.0)),
        periapsis_distance: h.q,
        time_to_periapsis,
        deflection: orbit::hyperbolic_deflection(mu, v_inf, b),
        eccentricity: h.e,
    })
}

/// Default start distance: well outside the system, and outside periapsis.
pub fn default_start_distance(system: &NBodySystem, mu: f64, v_inf: f64, b: f64) -> f64 {
    let extent = system.extent();
    let radius = system.bodies().iter().map(|b| b.radius()).fold(0.0, f64::max);
    let q = orbit::hyperbolic_periapsis(mu, v_inf, b);
    (10.0 * extent).max(100.0 * radius).max(3.0 * q).max(3.0 * b)
}

/// Create the projectile body on its incoming hyperbola and a tracker for it.
pub fn inject(system: &mut NBodySystem, spec: &ProjectileSpec) -> Result<(BodyId, ProjectileTracker), ProjectileError> {
    if system.is_empty() {
        return Err(ProjectileError::EmptySystem);
    }
    let m_proj = spec.structure.mass();
    let m_rest = system.total_mass();
    let mu = system.g * (m_rest + m_proj);
    let d = spec.start_distance.unwrap_or_else(|| default_start_distance(system, mu, spec.v_infinity, spec.impact_parameter));
    let (r, v) = hyperbolic_approach_state(mu, spec.v_infinity, spec.impact_parameter, d, spec.direction, spec.offset_direction)?;
    let com = system.center_of_mass();
    let vcom = system.com_velocity();
    let original_ids = system.ids();
    let body = Body::new(spec.name.clone(), spec.structure.clone(), com + r, vcom + v).with_kind(BodyKind::Projectile);
    let id = system.add(body);
    let mut tracker = ProjectileTracker::new(system, id, original_ids, d, spec.v_infinity, spec.direction.normalized());
    tracker.update(system, &[]);
    Ok((id, tracker))
}

/// Create a projectile at an explicit position and velocity (a launch from
/// inside the system rather than an arrival from infinity) and a tracker for
/// it. The tracker's "incoming" excess speed is the one implied by the launch
/// energy (zero when the launch is already bound), and an escape is declared
/// once the projectile has receded beyond the larger of its launch distance
/// and three times the system's extent.
pub fn inject_at(
    system: &mut NBodySystem,
    name: impl Into<String>,
    structure: SphericalStructure,
    position: Vec3,
    velocity: Vec3,
) -> Result<(BodyId, ProjectileTracker), ProjectileError> {
    if system.is_empty() {
        return Err(ProjectileError::EmptySystem);
    }
    let extent = system.extent();
    let original_ids = system.ids();
    let body = Body::new(name, structure, position, velocity).with_kind(BodyKind::Projectile);
    let id = system.add(body);
    let idx = system.index_of(id).expect("just added");
    let (eps, r, v, _) = specific_energy_relative(system, idx).ok_or(ProjectileError::EmptySystem)?;
    let v_inf = if eps > 0.0 { (2.0 * eps).sqrt() } else { 0.0 };
    let dir = v.try_normalized().unwrap_or(Vec3::X);
    let start_distance = r.norm().max(3.0 * extent);
    let mut tracker = ProjectileTracker::new(system, id, original_ids, start_distance, v_inf, dir);
    tracker.launched = true;
    // A launch that is not moving inward has, in effect, already had its periapsis.
    if r.dot(v) >= 0.0 {
        tracker.periapsis_passages = 1;
    }
    tracker.update(system, &[]);
    Ok((id, tracker))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Phase {
    /// Before the first periapsis passage.
    Approaching,
    /// Past periapsis, unbound and moving outward, but not yet past the start distance.
    Receding,
    /// Past periapsis, unbound, but falling back in (multi-body interplay).
    Interacting,
    /// Within a few Hill radii of a secondary body; verdicts are suspended.
    Encounter,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Outcome {
    InFlight { phase: Phase },
    Captured { orbits: u32, elements: OrbitalElements },
    Ejected { v_infinity_out: f64, deflection: f64, specific_energy_change: f64 },
    Merged { into: String, time: f64 },
    Undetermined,
}

impl Outcome {
    pub fn is_terminal(&self) -> bool {
        matches!(self, Outcome::Captured { .. } | Outcome::Ejected { .. } | Outcome::Merged { .. })
    }

    pub fn class(&self) -> OutcomeClass {
        match self {
            Outcome::Captured { .. } => OutcomeClass::Captured,
            Outcome::Ejected { .. } => OutcomeClass::Ejected,
            Outcome::Merged { .. } => OutcomeClass::Merged,
            Outcome::InFlight { .. } | Outcome::Undetermined => OutcomeClass::Undetermined,
        }
    }

    pub fn label(&self) -> String {
        match self {
            Outcome::InFlight { phase } => format!("in flight ({})", match phase {
                Phase::Approaching => "approaching",
                Phase::Receding => "receding",
                Phase::Interacting => "interacting",
                Phase::Encounter => "close encounter",
            }),
            Outcome::Captured { orbits, elements } => format!(
                "CAPTURED: bound orbit a={:.4e} m, e={:.3}, period={:.4e} s ({} full orbits so far)",
                elements.elements.a,
                elements.elements.e,
                elements.period.unwrap_or(f64::NAN),
                orbits
            ),
            Outcome::Ejected { v_infinity_out, deflection, specific_energy_change } => format!(
                "EJECTED: leaves with v_inf={:.4e} m/s, deflected by {:.2}°, Δε={:+.4e} J/kg",
                v_infinity_out,
                deflection.to_degrees(),
                specific_energy_change
            ),
            Outcome::Merged { into, time } => format!("MERGED into {into} at t={time:.4e} s"),
            Outcome::Undetermined => "undetermined".to_string(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OutcomeClass {
    Captured,
    Ejected,
    Merged,
    Undetermined,
}

impl OutcomeClass {
    pub fn symbol(self) -> char {
        match self {
            OutcomeClass::Captured => 'C',
            OutcomeClass::Ejected => '.',
            OutcomeClass::Merged => '#',
            OutcomeClass::Undetermined => '?',
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            OutcomeClass::Captured => "captured",
            OutcomeClass::Ejected => "ejected",
            OutcomeClass::Merged => "merged",
            OutcomeClass::Undetermined => "undetermined",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProjectileTracker {
    pub id: BodyId,
    pub name: String,
    pub original_ids: Vec<BodyId>,
    pub start_distance: f64,
    pub start_time: f64,
    pub incoming_v_infinity: f64,
    pub incoming_direction: Vec3,
    pub initial_specific_energy: f64,
    pub specific_energy: f64,
    pub distance: f64,
    pub closest_approach: f64,
    pub closest_approach_time: f64,
    /// Nearest any surface ever got (name, surface-to-surface gap in m).
    pub closest_body: Option<(String, f64)>,
    pub periapsis_passages: u32,
    pub apoapsis_passages: u32,
    pub outcome: Outcome,
    /// Original bodies absorbed by the projectile.
    pub absorbed: Vec<String>,
    /// Original bodies that were bound at injection and are now unbound.
    pub displaced: Vec<String>,
    /// The secondary whose Hill sphere the projectile is currently inside.
    pub encounter_with: Option<String>,
    /// `true` for a projectile launched from an explicit position inside the
    /// system ([`inject_at`]) rather than arriving from infinity.
    #[serde(default)]
    pub launched: bool,
    initially_bound: Vec<(BodyId, bool)>,
    last_radial_sign: i8,
}

/// Specific energy of body `idx` relative to all the other bodies.
fn specific_energy_relative(system: &NBodySystem, idx: usize) -> Option<(f64, Vec3, Vec3, f64)> {
    let bodies = system.bodies();
    let me = &bodies[idx];
    let mut m_rest = 0.0;
    let mut com = Vec3::ZERO;
    let mut vcom = Vec3::ZERO;
    let mut potential = 0.0;
    for (j, b) in bodies.iter().enumerate() {
        if j == idx {
            continue;
        }
        m_rest += b.mass();
        com += b.position * b.mass();
        vcom += b.velocity * b.mass();
        let d = me.position.distance(b.position);
        if d > 0.0 {
            potential -= system.g * b.mass() / d;
        }
    }
    if m_rest <= 0.0 {
        return None;
    }
    com /= m_rest;
    vcom /= m_rest;
    let r = me.position - com;
    let v = me.velocity - vcom;
    Some((0.5 * v.norm_sq() + potential, r, v, m_rest))
}

impl ProjectileTracker {
    fn new(system: &NBodySystem, id: BodyId, original_ids: Vec<BodyId>, start_distance: f64, v_inf: f64, direction: Vec3) -> Self {
        let initially_bound = original_ids
            .iter()
            .map(|&oid| {
                let idx = system.index_of(oid).expect("original id present");
                (oid, specific_energy_relative(system, idx).map(|e| e.0 < 0.0).unwrap_or(false))
            })
            .collect();
        let name = system.get(id).map(|b| b.name.clone()).unwrap_or_default();
        ProjectileTracker {
            id,
            name,
            original_ids,
            start_distance,
            start_time: system.time,
            incoming_v_infinity: v_inf,
            incoming_direction: direction,
            initial_specific_energy: 0.5 * v_inf * v_inf,
            specific_energy: 0.5 * v_inf * v_inf,
            distance: start_distance,
            closest_approach: f64::INFINITY,
            closest_approach_time: system.time,
            closest_body: None,
            periapsis_passages: 0,
            apoapsis_passages: 0,
            outcome: Outcome::InFlight { phase: Phase::Approaching },
            absorbed: Vec::new(),
            displaced: Vec::new(),
            encounter_with: None,
            launched: false,
            initially_bound,
            last_radial_sign: 0,
        }
    }

    /// Refresh the tracker from the current system state and any collision
    /// events that happened during the last step.
    pub fn update(&mut self, system: &NBodySystem, events: &[CollisionEvent]) -> &Outcome {
        if self.outcome.is_terminal() {
            return &self.outcome;
        }
        for e in events {
            if let CollisionOutcome::Merged { survivor, removed } = &e.outcome {
                if *removed == self.id {
                    let into = if e.a == *survivor { e.a_name.clone() } else { e.b_name.clone() };
                    self.outcome = Outcome::Merged { into, time: e.time };
                    return &self.outcome;
                }
                if *survivor == self.id {
                    let other = if e.a == *removed { e.a_name.clone() } else { e.b_name.clone() };
                    self.absorbed.push(other);
                }
            }
        }
        let Some(idx) = system.index_of(self.id) else {
            self.outcome = Outcome::Undetermined;
            return &self.outcome;
        };
        let Some((eps, r, v, m_rest)) = specific_energy_relative(system, idx) else {
            self.outcome = Outcome::Undetermined;
            return &self.outcome;
        };
        if self.last_radial_sign == 0 {
            self.initial_specific_energy = eps;
        }
        self.specific_energy = eps;
        let dist = r.norm();
        self.distance = dist;
        if dist < self.closest_approach {
            self.closest_approach = dist;
            self.closest_approach_time = system.time;
        }
        let me = &system.bodies()[idx];
        for (j, b) in system.bodies().iter().enumerate() {
            if j == idx {
                continue;
            }
            let gap = me.position.distance(b.position) - me.radius() - b.radius();
            if self.closest_body.as_ref().map(|(_, g)| gap < *g).unwrap_or(true) {
                self.closest_body = Some((b.name.clone(), gap));
            }
        }

        // Inside a few Hill radii of a secondary the projectile's energy
        // relative to the system is dominated by that body's local well and
        // its barycentric radial motion is meaningless, so verdicts and
        // periapsis/apoapsis counting are suspended until it is clear again.
        self.encounter_with = None;
        if let Some(primary) = system.bodies().iter().enumerate().filter(|(j, _)| *j != idx).max_by(|a, b| a.1.mass().partial_cmp(&b.1.mass()).unwrap_or(std::cmp::Ordering::Equal)) {
            let (pj, p) = primary;
            for (j, b) in system.bodies().iter().enumerate() {
                if j == idx || j == pj || b.mass() <= 0.0 {
                    continue;
                }
                let a_j = b.position.distance(p.position);
                let hill = a_j * (b.mass() / (3.0 * p.mass())).cbrt();
                if me.position.distance(b.position) < 3.0 * hill {
                    self.encounter_with = Some(b.name.clone());
                    break;
                }
            }
        }
        if self.encounter_with.is_some() {
            self.outcome = Outcome::InFlight { phase: Phase::Encounter };
            return &self.outcome;
        }

        let radial = r.dot(v);
        let sign: i8 = if radial > 0.0 { 1 } else if radial < 0.0 { -1 } else { 0 };
        if sign != 0 {
            if self.last_radial_sign < 0 && sign > 0 {
                self.periapsis_passages += 1;
            } else if self.last_radial_sign > 0 && sign < 0 {
                self.apoapsis_passages += 1;
            }
            self.last_radial_sign = sign;
        }

        if system.len() <= 256 {
            self.displaced.clear();
            for &(oid, was_bound) in &self.initially_bound {
                if !was_bound {
                    continue;
                }
                if let Some(j) = system.index_of(oid) {
                    if let Some((e, ..)) = specific_energy_relative(system, j) {
                        if e > 0.0 {
                            self.displaced.push(system.bodies()[j].name.clone());
                        }
                    }
                }
            }
        }

        // Energy is exchanged only during encounters, so once the projectile
        // is past periapsis, clear of every secondary, and has ε < 0, it is
        // bound: it will come back. Ejection is declared conservatively, once
        // it has receded beyond its start distance with ε > 0.
        self.outcome = if eps > 0.0 && radial > 0.0 && dist >= self.start_distance && self.periapsis_passages >= 1 {
            let v_out = (2.0 * eps).sqrt();
            let cos = v.normalized().dot(self.incoming_direction).clamp(-1.0, 1.0);
            Outcome::Ejected { v_infinity_out: v_out, deflection: cos.acos(), specific_energy_change: eps - self.initial_specific_energy }
        } else if eps < 0.0 && self.periapsis_passages >= 1 {
            let mu = system.g * (m_rest + me.mass());
            Outcome::Captured { orbits: self.apoapsis_passages, elements: orbit::elements_from_state(mu, r, v) }
        } else {
            let phase = if self.periapsis_passages == 0 {
                Phase::Approaching
            } else if radial > 0.0 {
                Phase::Receding
            } else {
                Phase::Interacting
            };
            Outcome::InFlight { phase }
        };
        &self.outcome
    }

    /// Osculating elements of the projectile about the rest of the system.
    pub fn current_elements(&self, system: &NBodySystem) -> Option<OrbitalElements> {
        let idx = system.index_of(self.id)?;
        let (_, r, v, m_rest) = specific_energy_relative(system, idx)?;
        let mu = system.g * (m_rest + system.bodies()[idx].mass());
        Some(orbit::elements_from_state(mu, r, v))
    }

    pub fn summary(&self) -> String {
        let mut s = format!("{}: {}", self.name, self.outcome.label());
        s.push_str(&format!(
            "\n  closest approach to barycentre {:.4e} m at t={:.4e} s; periapsis passages {}, apoapsis passages {}",
            self.closest_approach, self.closest_approach_time, self.periapsis_passages, self.apoapsis_passages
        ));
        if let Some((name, gap)) = &self.closest_body {
            s.push_str(&format!("\n  closest surface gap: {gap:.4e} m to {name}"));
        }
        s.push_str(&format!("\n  specific energy now {:+.4e} J/kg (started {:+.4e})", self.specific_energy, self.initial_specific_energy));
        if !self.absorbed.is_empty() {
            s.push_str(&format!("\n  absorbed: {}", self.absorbed.join(", ")));
        }
        if !self.displaced.is_empty() {
            s.push_str(&format!("\n  unbound by the encounter: {}", self.displaced.join(", ")));
        }
        s
    }
}

/// Outcome map over a grid of impact parameters and excess speeds.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SweepResult {
    pub impact_parameters: Vec<f64>,
    pub v_infinities: Vec<f64>,
    /// `grid[iv][ib]`.
    pub grid: Vec<Vec<OutcomeClass>>,
}

impl SweepResult {
    /// Fraction of the swept impact-parameter range (weighted by `2πb db`,
    /// i.e. by cross-section) that ends in each class, per speed.
    pub fn cross_section_fractions(&self) -> Vec<[f64; 4]> {
        self.grid
            .iter()
            .map(|row| {
                let mut w = [0.0; 4];
                let mut total = 0.0;
                for (ib, class) in row.iter().enumerate() {
                    let b = self.impact_parameters[ib];
                    let db = if self.impact_parameters.len() > 1 {
                        let lo = if ib > 0 { 0.5 * (b + self.impact_parameters[ib - 1]) } else { 0.0 };
                        let hi = if ib + 1 < self.impact_parameters.len() {
                            0.5 * (b + self.impact_parameters[ib + 1])
                        } else {
                            b + 0.5 * (b - self.impact_parameters[ib - 1])
                        };
                        hi * hi - lo * lo
                    } else {
                        1.0
                    };
                    let k = match class {
                        OutcomeClass::Captured => 0,
                        OutcomeClass::Ejected => 1,
                        OutcomeClass::Merged => 2,
                        OutcomeClass::Undetermined => 3,
                    };
                    w[k] += db;
                    total += db;
                }
                if total > 0.0 {
                    for x in &mut w {
                        *x /= total;
                    }
                }
                w
            })
            .collect()
    }

    pub fn ascii(&self) -> String {
        let mut s = String::new();
        s.push_str("rows: v_inf (m/s), columns: impact parameter (m); C=captured .=ejected #=merged ?=undetermined\n");
        for (iv, row) in self.grid.iter().enumerate() {
            s.push_str(&format!("{:>12.4e} | ", self.v_infinities[iv]));
            for c in row {
                s.push(c.symbol());
            }
            s.push('\n');
        }
        s.push_str(&format!(
            "{:>12} | b from {:.4e} to {:.4e}\n",
            "",
            self.impact_parameters.first().copied().unwrap_or(0.0),
            self.impact_parameters.last().copied().unwrap_or(0.0)
        ));
        s
    }
}

/// Run one projectile experiment per `(v_inf, b)` pair on a copy of `base`.
pub fn sweep(
    base: &NBodySystem,
    config: &SimulationConfig,
    template: &ProjectileSpec,
    impact_parameters: &[f64],
    v_infinities: &[f64],
    max_time: f64,
    mut progress: Option<&mut dyn FnMut(usize, usize)>,
) -> SweepResult {
    let total = impact_parameters.len() * v_infinities.len();
    let mut done = 0;
    let mut grid = Vec::with_capacity(v_infinities.len());
    for &v in v_infinities {
        let mut row = Vec::with_capacity(impact_parameters.len());
        for &b in impact_parameters {
            let mut sim = Simulation::new(base.clone(), config.clone());
            let mut spec = template.clone();
            spec.v_infinity = v;
            spec.impact_parameter = b;
            let class = match sim.inject_projectile(&spec) {
                Ok(_) => sim.run_projectile(max_time).class(),
                Err(_) => OutcomeClass::Undetermined,
            };
            row.push(class);
            done += 1;
            if let Some(p) = progress.as_deref_mut() {
                p(done, total);
            }
        }
        grid.push(row);
    }
    SweepResult { impact_parameters: impact_parameters.to_vec(), v_infinities: v_infinities.to_vec(), grid }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collision::CollisionMode;
    use crate::force::ForceMethod;
    use crate::integrator::IntegratorKind;
    use crate::material::Material;
    use crate::simulation::TimestepPolicy;
    use approx::assert_relative_eq;

    #[test]
    fn approach_state_has_correct_energy_and_angular_momentum() {
        let (mu, v_inf, b, d) = (1.0, 0.7, 2.5, 40.0);
        let dir = Vec3::new(1.0, 0.3, -0.2);
        let (r, v) = hyperbolic_approach_state(mu, v_inf, b, d, dir, None).unwrap();
        assert_relative_eq!(r.norm(), d, max_relative = 1e-12);
        let energy = 0.5 * v.norm_sq() - mu / r.norm();
        assert_relative_eq!(energy, 0.5 * v_inf * v_inf, max_relative = 1e-12);
        assert_relative_eq!(r.cross(v).norm(), b * v_inf, max_relative = 1e-12);
        // Still approaching, and roughly along the requested direction.
        assert!(r.dot(v) < 0.0);
        assert!(v.normalized().dot(dir.normalized()) > 0.95);
        // Angular momentum perpendicular to both the direction and the offset.
        let off = Vec3::Z.cross(dir.normalized()).normalized();
        let h = r.cross(v).normalized();
        assert!(h.dot(dir.normalized()).abs() < 1e-12);
        assert!(h.dot(off).abs() < 1e-12);
        // The undeflected path is displaced towards +offset.
        let perp = r.reject_from_unit(dir.normalized());
        assert!(perp.dot(off) > 0.0);
    }

    #[test]
    fn approach_geometry_predicts_periapsis_time_and_point() {
        // Integrate the two-body problem and compare with the analytic prediction.
        let (mu, v_inf, b, d) = (1.0, 0.8, 1.7, 30.0);
        let geo = approach_geometry(mu, v_inf, b, d, Vec3::new(0.2, 1.0, 0.0), None).unwrap();
        assert_relative_eq!(geo.periapsis_distance, orbit::hyperbolic_periapsis(mu, v_inf, b), max_relative = 1e-12);
        assert_relative_eq!(geo.periapsis.norm(), geo.periapsis_distance, max_relative = 1e-12);
        let mut s = crate::integrator::State::new(vec![geo.position], vec![geo.velocity]);
        let mut integ = crate::integrator::IntegratorKind::Pefrl.build();
        let dt = 1e-4;
        let (mut t, mut best, mut best_t, mut best_r) = (0.0, f64::INFINITY, 0.0, Vec3::ZERO);
        let mut accel = |p: &[Vec3], out: &mut Vec<Vec3>| {
            out.clear();
            out.push(-p[0] * (mu / p[0].norm().powi(3)));
        };
        while t < 2.0 * geo.time_to_periapsis {
            integ.step(&mut s, dt, &mut accel);
            t += dt;
            let r = s.pos[0].norm();
            if r < best {
                best = r;
                best_t = t;
                best_r = s.pos[0];
            }
        }
        assert_relative_eq!(best_t, geo.time_to_periapsis, max_relative = 1e-3);
        assert!((best_r - geo.periapsis).norm() < 1e-3 * geo.periapsis_distance);
        // Radial case: pure infall time.
        let radial = approach_geometry(mu, v_inf, 0.0, d, Vec3::X, None).unwrap();
        assert!(radial.time_to_periapsis > 0.0 && radial.time_to_periapsis < d / v_inf);
    }

    #[test]
    fn radial_approach_state() {
        let (r, v) = hyperbolic_approach_state(1.0, 1.0, 0.0, 10.0, Vec3::X, None).unwrap();
        assert_eq!(r, Vec3::new(-10.0, 0.0, 0.0));
        assert_relative_eq!(v.x, (1.0f64 + 0.2).sqrt(), max_relative = 1e-12);
    }

    #[test]
    fn single_body_flyby_is_ejected_with_rutherford_deflection() {
        let mut s = NBodySystem::new(1.0);
        s.add(Body::uniform("M", Material::Iron, 1.0, Vec3::ZERO, Vec3::ZERO).unwrap());
        let proj = SphericalStructure::uniform_with_mass(Material::Rock, 1e-9).unwrap();
        let (v_inf, b) = (1.0, 2.0);
        let spec = ProjectileSpec::new("p", proj, v_inf, b).with_start_distance(60.0);
        let cfg = SimulationConfig {
            integrator: IntegratorKind::Pefrl,
            timestep: TimestepPolicy::Adaptive { eta: 0.01, dt_min: 1e-6, dt_max: 0.05 },
            force: ForceMethod::Direct,
            collisions: CollisionMode::Physical,
        };
        let mut sim = Simulation::new(s, cfg);
        sim.inject_projectile(&spec).unwrap();
        let outcome = sim.run_projectile(1e4);
        match outcome {
            Outcome::Ejected { v_infinity_out, deflection, .. } => {
                assert_relative_eq!(v_infinity_out, v_inf, max_relative = 1e-4);
                let expected = orbit::hyperbolic_deflection(1.0 + 1e-9, v_inf, b);
                assert_relative_eq!(deflection, expected, max_relative = 2e-3);
            }
            other => panic!("expected ejection, got {other:?}"),
        }
        let t = sim.trackers.last().unwrap();
        let q = orbit::hyperbolic_periapsis(1.0, v_inf, b);
        assert_relative_eq!(t.closest_approach, q, max_relative = 1e-3);
        assert!(t.displaced.is_empty());
    }

    #[test]
    fn launched_projectiles_are_classified() {
        let cfg = SimulationConfig {
            integrator: IntegratorKind::Pefrl,
            timestep: TimestepPolicy::Adaptive { eta: 0.01, dt_min: 1e-6, dt_max: 0.05 },
            force: ForceMethod::Direct,
            collisions: CollisionMode::Physical,
        };
        let proj = SphericalStructure::uniform_with_mass(Material::Rock, 1e-9).unwrap();
        let mut s = NBodySystem::new(1.0);
        s.add(Body::uniform("M", Material::Iron, 1.0, Vec3::ZERO, Vec3::ZERO).unwrap());
        // Below the escape speed: a bound orbit, reported with its elements.
        let mut sim = Simulation::new(s.clone(), cfg.clone());
        let id = sim.inject_projectile_at("slow", proj.clone(), Vec3::new(2.0, 0.0, 0.0), Vec3::new(0.0, 0.5, 0.0)).unwrap();
        assert!(sim.tracker(id).unwrap().launched);
        let outcome = sim.run_projectile(200.0);
        match outcome {
            Outcome::Captured { elements, .. } => {
                let expected = elements_from_launch(1.0, Vec3::new(2.0, 0.0, 0.0), Vec3::new(0.0, 0.5, 0.0));
                assert_relative_eq!(elements.elements.a, expected, max_relative = 1e-3);
            }
            other => panic!("expected a bound orbit, got {other:?}"),
        }
        // Above it: escapes, and the tracker reports the excess speed.
        let mut sim = Simulation::new(s, cfg);
        sim.inject_projectile_at("fast", proj, Vec3::new(2.0, 0.0, 0.0), Vec3::new(0.0, 1.5, 0.0)).unwrap();
        let outcome = sim.run_projectile(500.0);
        match outcome {
            Outcome::Ejected { v_infinity_out, .. } => {
                let v_inf = (1.5f64 * 1.5 - 2.0 / 2.0).sqrt();
                assert_relative_eq!(v_infinity_out, v_inf, max_relative = 1e-3);
            }
            other => panic!("expected an escape, got {other:?}"),
        }
    }

    fn elements_from_launch(mu: f64, r: Vec3, v: Vec3) -> f64 {
        orbit::elements_from_state(mu, r, v).elements.a
    }

    #[test]
    fn head_on_projectile_merges() {
        let mut s = NBodySystem::new(1.0);
        s.add(Body::uniform("M", Material::Iron, 1.0, Vec3::ZERO, Vec3::ZERO).unwrap());
        let proj = SphericalStructure::uniform_with_mass(Material::Rock, 1e-3).unwrap();
        let spec = ProjectileSpec::new("p", proj, 0.5, 0.0).with_start_distance(20.0);
        let cfg = SimulationConfig {
            integrator: IntegratorKind::Leapfrog,
            timestep: TimestepPolicy::Adaptive { eta: 0.01, dt_min: 1e-7, dt_max: 0.05 },
            force: ForceMethod::Direct,
            collisions: CollisionMode::Physical,
        };
        let mut sim = Simulation::new(s, cfg);
        sim.inject_projectile(&spec).unwrap();
        let outcome = sim.run_projectile(1e3);
        assert!(matches!(outcome, Outcome::Merged { .. }), "{outcome:?}");
        assert_eq!(sim.system.len(), 1);
        assert_relative_eq!(sim.system.total_mass(), 1.0 + 1e-3, max_relative = 1e-12);
    }
}
