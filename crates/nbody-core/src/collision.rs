//! Finite-size collisions between spherical bodies.
//!
//! Two bodies collide when their surfaces overlap, `|r_ij| < R_i + R_j`.
//! Three resolutions are available:
//!
//! * **Merge**: perfectly inelastic. Mass and momentum are conserved, the
//!   survivor sits at the centre of mass and the internal structure is the
//!   differentiated union of both bodies (see
//!   [`SphericalStructure::merged_with`](crate::structure::SphericalStructure::merged_with)).
//! * **Bounce**: an impulsive normal collision with the geometric-mean
//!   coefficient of restitution of the two materials, plus a positional
//!   correction that removes the overlap. Momentum is conserved exactly;
//!   kinetic energy along the normal is reduced by `1 − e²`.
//! * **Physical** (default): the regime map of Leinhardt & Stewart (2012,
//!   *ApJ* 745, 79) without fragmentation. If the impact speed is below the
//!   mutual escape speed the bodies cannot separate again, so they merge;
//!   faster impacts merge when they are not grazing and *hit-and-run*
//!   (bounce) when they are. Catastrophic disruption is not modelled: a fast
//!   head-on impact accretes instead of shattering.

use crate::body::{Body, BodyId, BodyKind};
use crate::system::NBodySystem;
use crate::vector::Vec3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum CollisionMode {
    /// Bodies pass through each other (pure point masses).
    None,
    Merge,
    Bounce,
    #[default]
    Physical,
}

impl CollisionMode {
    pub const ALL: [CollisionMode; 4] = [CollisionMode::None, CollisionMode::Merge, CollisionMode::Bounce, CollisionMode::Physical];

    pub fn name(self) -> &'static str {
        match self {
            CollisionMode::None => "none",
            CollisionMode::Merge => "merge",
            CollisionMode::Bounce => "bounce",
            CollisionMode::Physical => "physical",
        }
    }

    pub fn parse(s: &str) -> Option<CollisionMode> {
        Some(match s.trim().to_ascii_lowercase().as_str() {
            "none" | "off" | "ignore" => CollisionMode::None,
            "merge" | "accrete" => CollisionMode::Merge,
            "bounce" | "elastic" => CollisionMode::Bounce,
            "physical" | "auto" => CollisionMode::Physical,
            _ => return None,
        })
    }
}

impl std::fmt::Display for CollisionMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CollisionOutcome {
    Merged { survivor: BodyId, removed: BodyId },
    Bounced { restitution: f64 },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CollisionEvent {
    pub time: f64,
    pub a: BodyId,
    pub a_name: String,
    pub b: BodyId,
    pub b_name: String,
    /// Relative speed at contact, m/s.
    pub impact_speed: f64,
    /// Mutual escape speed `sqrt(2G(m_a+m_b)/(R_a+R_b))`, m/s.
    pub escape_speed: f64,
    /// Impact angle: 0 is head-on, π/2 is a graze, rad.
    pub impact_angle: f64,
    /// Change of the system's mechanical energy `K + U` caused by the event, J.
    ///
    /// A bounce only removes kinetic energy. A merger also deletes the pair's
    /// mutual potential term `−G m_a m_b / d` (it becomes internal binding
    /// energy of the remnant), so its `energy_change` can have either sign;
    /// what is always dissipated is the kinetic energy of relative motion
    /// `½ μ v_rel²`.
    pub energy_change: f64,
    pub outcome: CollisionOutcome,
}

impl CollisionEvent {
    pub fn describe(&self) -> String {
        match &self.outcome {
            CollisionOutcome::Merged { survivor, .. } => {
                let keep = if *survivor == self.a { &self.a_name } else { &self.b_name };
                format!(
                    "{} and {} merged into {} (impact {:.3} km/s vs escape {:.3} km/s, angle {:.0}°)",
                    self.a_name,
                    self.b_name,
                    keep,
                    self.impact_speed / 1e3,
                    self.escape_speed / 1e3,
                    self.impact_angle.to_degrees()
                )
            }
            CollisionOutcome::Bounced { restitution } => format!(
                "{} and {} bounced, e={:.2} (impact {:.3} km/s vs escape {:.3} km/s, angle {:.0}°)",
                self.a_name,
                self.b_name,
                restitution,
                self.impact_speed / 1e3,
                self.escape_speed / 1e3,
                self.impact_angle.to_degrees()
            ),
        }
    }
}

/// Indices of all overlapping pairs.
pub fn detect_overlaps(bodies: &[Body]) -> Vec<(usize, usize)> {
    let n = bodies.len();
    let mut out = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let sum = bodies[i].radius() + bodies[j].radius();
            if bodies[i].position.distance_sq(bodies[j].position) < sum * sum {
                out.push((i, j));
            }
        }
    }
    out
}

struct Contact {
    normal: Vec3,
    impact_speed: f64,
    escape_speed: f64,
    impact_angle: f64,
    approaching: bool,
}

fn contact(a: &Body, b: &Body, g: f64) -> Contact {
    let d = b.position - a.position;
    let normal = d.try_normalized().unwrap_or_else(|| (b.velocity - a.velocity).try_normalized().unwrap_or(Vec3::X));
    let v_rel = b.velocity - a.velocity;
    let vn = v_rel.dot(normal);
    let speed = v_rel.norm();
    let impact_angle = if speed > 0.0 { (v_rel.cross(normal).norm() / speed).clamp(0.0, 1.0).asin() } else { 0.0 };
    Contact {
        normal,
        impact_speed: speed,
        escape_speed: (2.0 * g * (a.mass() + b.mass()) / (a.radius() + b.radius())).sqrt(),
        impact_angle,
        approaching: vn < 0.0,
    }
}

/// Decide merge vs bounce for a contact in `Physical` mode.
fn physical_merges(a: &Body, b: &Body, c: &Contact) -> bool {
    if c.impact_speed <= c.escape_speed {
        return true;
    }
    // Leinhardt & Stewart (2012): grazing when sin(θ) > R_target / (R_target + R_projectile).
    let (rt, rp) = if a.radius() >= b.radius() { (a.radius(), b.radius()) } else { (b.radius(), a.radius()) };
    let b_crit = rt / (rt + rp);
    c.impact_angle.sin() <= b_crit
}

fn merge(system: &mut NBodySystem, i: usize, j: usize) -> (BodyId, BodyId) {
    let (ia, ib) = if system.bodies()[i].mass() >= system.bodies()[j].mass() { (i, j) } else { (j, i) };
    let removed = system.remove_at(ib);
    // Removing shifts indices above `ib` down by one.
    let ia = if ib < ia { ia - 1 } else { ia };
    let survivor = &mut system.bodies_mut()[ia];
    let m_total = survivor.mass() + removed.mass();
    let com = (survivor.position * survivor.mass() + removed.position * removed.mass()) / m_total;
    let vel = (survivor.velocity * survivor.mass() + removed.velocity * removed.mass()) / m_total;
    let structure = survivor.structure().merged_with(removed.structure());
    survivor.set_structure(structure);
    survivor.position = com;
    survivor.velocity = vel;
    if survivor.kind != BodyKind::Star && removed.kind != BodyKind::Star {
        if matches!(survivor.kind, BodyKind::Projectile | BodyKind::Asteroid | BodyKind::Comet | BodyKind::Particle) {
            survivor.kind = BodyKind::Remnant;
        }
    } else if removed.kind == BodyKind::Star {
        survivor.kind = BodyKind::Star;
    }
    (survivor.id(), removed.id())
}

fn bounce(system: &mut NBodySystem, i: usize, j: usize, c: &Contact, restitution: f64) {
    let bodies = system.bodies_mut();
    let (ma, mb) = (bodies[i].mass(), bodies[j].mass());
    let inv = 1.0 / ma + 1.0 / mb;
    let v_rel = bodies[j].velocity - bodies[i].velocity;
    let vn = v_rel.dot(c.normal);
    // Impulse magnitude along the normal for restitution e.
    let impulse = -(1.0 + restitution) * vn / inv;
    bodies[i].velocity -= c.normal * (impulse / ma);
    bodies[j].velocity += c.normal * (impulse / mb);
    // Positional correction: push apart along the normal by the overlap, split by inverse mass.
    let overlap = bodies[i].radius() + bodies[j].radius() - bodies[i].position.distance(bodies[j].position);
    if overlap > 0.0 {
        let total = inv;
        bodies[i].position -= c.normal * (overlap * (1.0 / ma) / total);
        bodies[j].position += c.normal * (overlap * (1.0 / mb) / total);
    }
}

/// Detect and resolve all overlaps according to `mode`. Returns the events.
///
/// Energy bookkeeping is done per event using the full system energy, so
/// callers can separate physical dissipation from integration error.
pub fn resolve(system: &mut NBodySystem, mode: CollisionMode) -> Vec<CollisionEvent> {
    let mut events = Vec::new();
    if mode == CollisionMode::None {
        return events;
    }
    // Repeat until no approaching overlaps remain (a merge can create new overlaps).
    for _pass in 0..16 {
        let pairs = detect_overlaps(system.bodies());
        let mut acted = false;
        for (i, j) in pairs {
            if i >= system.len() || j >= system.len() {
                break; // indices invalidated by a merge; re-detect
            }
            let (a, b) = (&system.bodies()[i], &system.bodies()[j]);
            let c = contact(a, b, system.g);
            let merge_it = match mode {
                CollisionMode::Merge => true,
                CollisionMode::Bounce => false,
                CollisionMode::Physical => physical_merges(a, b, &c),
                CollisionMode::None => unreachable!(),
            };
            if !merge_it && !c.approaching {
                continue; // overlapping but already separating
            }
            let (a_id, a_name, b_id, b_name) = (a.id(), a.name.clone(), b.id(), b.name.clone());
            let e_before = system.total_energy();
            let outcome = if merge_it {
                let (survivor, removed) = merge(system, i, j);
                CollisionOutcome::Merged { survivor, removed }
            } else {
                let restitution = (system.bodies()[i].restitution() * system.bodies()[j].restitution()).sqrt();
                bounce(system, i, j, &c, restitution);
                CollisionOutcome::Bounced { restitution }
            };
            let e_after = system.total_energy();
            events.push(CollisionEvent {
                time: system.time,
                a: a_id,
                a_name,
                b: b_id,
                b_name,
                impact_speed: c.impact_speed,
                escape_speed: c.escape_speed,
                impact_angle: c.impact_angle,
                energy_change: e_after - e_before,
                outcome,
            });
            acted = true;
            if merge_it {
                break; // indices changed; re-detect from scratch
            }
        }
        if !acted {
            break;
        }
    }
    events
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::Material;
    use approx::assert_relative_eq;

    fn pair(mode_speed: f64) -> NBodySystem {
        let mut s = NBodySystem::new(1.0);
        s.add(Body::uniform("A", Material::Rock, 2.0, Vec3::ZERO, Vec3::ZERO).unwrap());
        let r = s.bodies()[0].radius();
        s.add(Body::uniform("B", Material::Ice, 1.0, Vec3::new(1.5 * r, 0.0, 0.0), Vec3::new(-mode_speed, 0.1, 0.0)).unwrap());
        s
    }

    #[test]
    fn merge_conserves_mass_and_momentum() {
        let mut s = pair(1.0);
        let p0 = s.momentum();
        let m0 = s.total_mass();
        let (a, b) = (&s.bodies()[0], &s.bodies()[1]);
        let d = a.position.distance(b.position);
        let mu_red = a.mass() * b.mass() / (a.mass() + b.mass());
        // Lost relative kinetic energy, plus the pair's potential term that
        // becomes internal binding energy of the remnant.
        let expected_de = -0.5 * mu_red * (b.velocity - a.velocity).norm_sq() + s.g * a.mass() * b.mass() / d;
        let ev = resolve(&mut s, CollisionMode::Merge);
        assert_eq!(ev.len(), 1);
        assert_eq!(s.len(), 1);
        assert_relative_eq!(s.total_mass(), m0, max_relative = 1e-12);
        assert!((s.momentum() - p0).norm() < 1e-12);
        assert_eq!(s.bodies()[0].name, "A");
        assert_eq!(s.bodies()[0].structure().layers().len(), 2);
        assert_relative_eq!(ev[0].energy_change, expected_de, max_relative = 1e-9);
    }

    #[test]
    fn bounce_conserves_momentum_and_dissipates_by_restitution() {
        let mut s = pair(1.0);
        let p0 = s.momentum();
        let vn0 = (s.bodies()[1].velocity - s.bodies()[0].velocity).x;
        let ev = resolve(&mut s, CollisionMode::Bounce);
        assert_eq!(ev.len(), 1);
        assert_eq!(s.len(), 2);
        assert!((s.momentum() - p0).norm() < 1e-12);
        let vn1 = (s.bodies()[1].velocity - s.bodies()[0].velocity).x;
        let e = (Material::Rock.restitution() * Material::Ice.restitution()).sqrt();
        assert_relative_eq!(vn1, -e * vn0, max_relative = 1e-12);
        // Separated after the positional correction.
        let d = s.bodies()[0].position.distance(s.bodies()[1].position);
        assert!(d >= s.bodies()[0].radius() + s.bodies()[1].radius() - 1e-9);
        // Separating bodies are left alone on a second pass.
        assert!(resolve(&mut s, CollisionMode::Bounce).is_empty());
    }

    #[test]
    fn physical_mode_merges_slow_and_bounces_fast_grazing() {
        let mut slow = pair(1e-6);
        let ev = resolve(&mut slow, CollisionMode::Physical);
        assert!(matches!(ev[0].outcome, CollisionOutcome::Merged { .. }));

        // Fast and grazing: approach mostly tangential.
        let mut s = NBodySystem::new(1.0);
        s.add(Body::uniform("A", Material::Iron, 2.0, Vec3::ZERO, Vec3::ZERO).unwrap());
        let r = s.bodies()[0].radius();
        let rb = Body::uniform("B", Material::Iron, 1.0, Vec3::ZERO, Vec3::ZERO).unwrap().radius();
        let sum = r + rb;
        s.add(
            Body::uniform("B", Material::Iron, 1.0, Vec3::new(0.0, 0.999 * sum, 0.0), Vec3::new(-100.0, -0.1, 0.0)).unwrap(),
        );
        let ev = resolve(&mut s, CollisionMode::Physical);
        assert_eq!(ev.len(), 1);
        assert!(matches!(ev[0].outcome, CollisionOutcome::Bounced { .. }), "{:?}", ev[0]);
        assert!(ev[0].impact_angle > 1.4);
    }

    #[test]
    fn none_mode_ignores_overlaps() {
        let mut s = pair(1.0);
        assert!(resolve(&mut s, CollisionMode::None).is_empty());
        assert_eq!(s.len(), 2);
    }
}
