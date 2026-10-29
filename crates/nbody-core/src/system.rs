//! The collection of bodies plus the constants that define the dynamics.

use crate::body::{Body, BodyId};
use crate::force;
use crate::orbit::{self, Elements};
use crate::vector::Vec3;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SystemError {
    #[error("no body named '{0}'")]
    NoSuchName(String),
    #[error("no body with id {0}")]
    NoSuchId(BodyId),
    #[error("a body cannot orbit itself")]
    SelfOrbit,
    #[error("invalid orbit: {0}")]
    BadOrbit(String),
}

/// All bodies, the gravitational constant, the softening length and the clock.
///
/// Every mutation through this API bumps a `generation` counter so that
/// callers caching derived data (the integrator's acceleration cache, GUI
/// trails) can detect edits cheaply.
#[derive(Debug, Clone)]
pub struct NBodySystem {
    bodies: Vec<Body>,
    /// Gravitational constant in the units of the state (SI by default).
    pub g: f64,
    /// Plummer softening length ε: pair force uses `(r² + ε²)^{3/2}`.
    pub softening: f64,
    /// Simulation clock.
    pub time: f64,
    next_id: u64,
    generation: u64,
}

impl NBodySystem {
    pub fn new(g: f64) -> Self {
        NBodySystem { bodies: Vec::new(), g, softening: 0.0, time: 0.0, next_id: 1, generation: 0 }
    }

    /// SI units with the CODATA value of `G`.
    pub fn si() -> Self {
        Self::new(crate::units::G_SI)
    }

    pub fn generation(&self) -> u64 {
        self.generation
    }

    /// Mark the system as modified (invalidates caches held by observers).
    pub fn touch(&mut self) {
        self.generation += 1;
    }

    /// Add a body; it receives a fresh id which is returned.
    pub fn add(&mut self, mut body: Body) -> BodyId {
        let id = BodyId(self.next_id);
        self.next_id += 1;
        body.set_id(id);
        self.bodies.push(body);
        self.touch();
        id
    }

    pub fn remove(&mut self, id: BodyId) -> Option<Body> {
        let idx = self.index_of(id)?;
        Some(self.remove_at(idx))
    }

    pub fn remove_at(&mut self, index: usize) -> Body {
        self.touch();
        self.bodies.remove(index)
    }

    pub fn clear(&mut self) {
        self.bodies.clear();
        self.touch();
    }

    pub fn len(&self) -> usize {
        self.bodies.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bodies.is_empty()
    }

    pub fn bodies(&self) -> &[Body] {
        &self.bodies
    }

    pub fn bodies_mut(&mut self) -> &mut [Body] {
        self.touch();
        &mut self.bodies
    }

    pub fn index_of(&self, id: BodyId) -> Option<usize> {
        self.bodies.iter().position(|b| b.id() == id)
    }

    pub fn get(&self, id: BodyId) -> Option<&Body> {
        self.bodies.iter().find(|b| b.id() == id)
    }

    pub fn get_mut(&mut self, id: BodyId) -> Option<&mut Body> {
        self.touch();
        self.bodies.iter_mut().find(|b| b.id() == id)
    }

    pub fn find(&self, name: &str) -> Option<&Body> {
        self.bodies.iter().find(|b| b.name.eq_ignore_ascii_case(name))
    }

    pub fn find_mut(&mut self, name: &str) -> Option<&mut Body> {
        self.touch();
        self.bodies.iter_mut().find(|b| b.name.eq_ignore_ascii_case(name))
    }

    pub fn id_of(&self, name: &str) -> Result<BodyId, SystemError> {
        self.find(name).map(|b| b.id()).ok_or_else(|| SystemError::NoSuchName(name.to_string()))
    }

    pub fn ids(&self) -> Vec<BodyId> {
        self.bodies.iter().map(|b| b.id()).collect()
    }

    /// The most massive body (the natural "primary").
    pub fn most_massive(&self) -> Option<&Body> {
        self.bodies.iter().max_by(|a, b| a.mass().partial_cmp(&b.mass()).unwrap_or(std::cmp::Ordering::Equal))
    }

    // ---- bulk state access -------------------------------------------------

    pub fn masses(&self) -> Vec<f64> {
        self.bodies.iter().map(|b| b.mass()).collect()
    }

    pub fn positions(&self) -> Vec<Vec3> {
        self.bodies.iter().map(|b| b.position).collect()
    }

    pub fn velocities(&self) -> Vec<Vec3> {
        self.bodies.iter().map(|b| b.velocity).collect()
    }

    /// Copy positions, velocities and masses into reusable buffers.
    pub fn fill_state(&self, pos: &mut Vec<Vec3>, vel: &mut Vec<Vec3>, mass: &mut Vec<f64>) {
        pos.clear();
        vel.clear();
        mass.clear();
        for b in &self.bodies {
            pos.push(b.position);
            vel.push(b.velocity);
            mass.push(b.mass());
        }
    }

    /// Write back integrated positions and velocities. Does not bump the
    /// generation: the integrator is the one observer that already knows.
    pub fn write_state(&mut self, pos: &[Vec3], vel: &[Vec3]) {
        debug_assert_eq!(pos.len(), self.bodies.len());
        for (b, (p, v)) in self.bodies.iter_mut().zip(pos.iter().zip(vel)) {
            b.position = *p;
            b.velocity = *v;
        }
    }

    // ---- integrals and bulk properties ------------------------------------

    pub fn total_mass(&self) -> f64 {
        self.bodies.iter().map(|b| b.mass()).sum()
    }

    pub fn center_of_mass(&self) -> Vec3 {
        let m = self.total_mass();
        if m <= 0.0 {
            return Vec3::ZERO;
        }
        self.bodies.iter().map(|b| b.position * b.mass()).sum::<Vec3>() / m
    }

    pub fn com_velocity(&self) -> Vec3 {
        let m = self.total_mass();
        if m <= 0.0 {
            return Vec3::ZERO;
        }
        self.bodies.iter().map(|b| b.velocity * b.mass()).sum::<Vec3>() / m
    }

    /// Shift to the barycentric frame (centre of mass at rest at the origin).
    pub fn move_to_barycentre(&mut self) {
        let r = self.center_of_mass();
        let v = self.com_velocity();
        for b in &mut self.bodies {
            b.position -= r;
            b.velocity -= v;
        }
        self.touch();
    }

    pub fn kinetic_energy(&self) -> f64 {
        self.bodies.iter().map(|b| b.kinetic_energy()).sum()
    }

    /// Plummer-softened potential energy, consistent with the force law.
    pub fn potential_energy(&self) -> f64 {
        let pos = self.positions();
        let mass = self.masses();
        force::potential_energy(&pos, &mass, self.g, self.softening)
    }

    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }

    pub fn momentum(&self) -> Vec3 {
        self.bodies.iter().map(|b| b.momentum()).sum()
    }

    /// Total angular momentum about the origin.
    pub fn angular_momentum(&self) -> Vec3 {
        self.bodies.iter().map(|b| b.angular_momentum(Vec3::ZERO)).sum()
    }

    /// Largest distance of any body from the centre of mass.
    pub fn extent(&self) -> f64 {
        let c = self.center_of_mass();
        self.bodies.iter().map(|b| b.position.distance(c)).fold(0.0, f64::max)
    }

    /// `G·M` of one body.
    pub fn mu_of(&self, id: BodyId) -> Result<f64, SystemError> {
        self.get(id).map(|b| self.g * b.mass()).ok_or(SystemError::NoSuchId(id))
    }

    // ---- orbit placement ---------------------------------------------------

    /// Put `id` on a circular orbit around `around` at its current distance.
    /// `normal` is the orbital angular-momentum direction (default +z).
    /// The two-body reduced problem is used, so `μ = G(M + m)`.
    pub fn set_circular_orbit(&mut self, id: BodyId, around: BodyId, normal: Vec3) -> Result<(), SystemError> {
        if id == around {
            return Err(SystemError::SelfOrbit);
        }
        let primary = self.get(around).ok_or(SystemError::NoSuchId(around))?.clone();
        let body = self.get(id).ok_or(SystemError::NoSuchId(id))?;
        let rel = body.position - primary.position;
        let r = rel.norm();
        if r <= 0.0 {
            return Err(SystemError::BadOrbit("bodies are at the same position".into()));
        }
        let n = normal.try_normalized().unwrap_or(Vec3::Z);
        let mut tangent = n.cross(rel);
        if tangent.norm() == 0.0 {
            tangent = rel.any_orthonormal();
        }
        let mu = self.g * (primary.mass() + body.mass());
        let v = (mu / r).sqrt();
        let vel = primary.velocity + tangent.normalized() * v;
        self.get_mut(id).expect("checked").velocity = vel;
        Ok(())
    }

    /// Place `id` on a Keplerian orbit (elements in SI / radians) about `around`.
    pub fn set_orbit(&mut self, id: BodyId, around: BodyId, elements: &Elements) -> Result<(), SystemError> {
        if id == around {
            return Err(SystemError::SelfOrbit);
        }
        let primary = self.get(around).ok_or(SystemError::NoSuchId(around))?.clone();
        let body = self.get(id).ok_or(SystemError::NoSuchId(id))?;
        let mu = self.g * (primary.mass() + body.mass());
        let (r, v) = orbit::state_from_elements(mu, elements).map_err(|e| SystemError::BadOrbit(e.to_string()))?;
        let b = self.get_mut(id).expect("checked");
        b.position = primary.position + r;
        b.velocity = primary.velocity + v;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::Material;
    use approx::assert_relative_eq;

    #[test]
    fn circular_orbit_velocity() {
        let mut s = NBodySystem::new(1.0);
        let a = s.add(Body::uniform("A", Material::Rock, 1.0, Vec3::ZERO, Vec3::ZERO).unwrap());
        let b = s.add(Body::uniform("B", Material::Rock, 1e-3, Vec3::new(2.0, 0.0, 0.0), Vec3::ZERO).unwrap());
        s.set_circular_orbit(b, a, Vec3::Z).unwrap();
        let v = s.get(b).unwrap().velocity;
        assert_relative_eq!(v.y, (1.001f64 / 2.0).sqrt(), max_relative = 1e-12);
        assert_eq!(v.x, 0.0);
    }

    #[test]
    fn barycentre_shift_zeroes_momentum() {
        let mut s = NBodySystem::new(1.0);
        s.add(Body::uniform("A", Material::Rock, 3.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap());
        s.add(Body::uniform("B", Material::Rock, 1.0, Vec3::new(-1.0, 2.0, 0.0), Vec3::new(1.0, 0.0, 0.0)).unwrap());
        s.move_to_barycentre();
        assert!(s.momentum().norm() < 1e-12);
        assert!(s.center_of_mass().norm() < 1e-12);
    }

    #[test]
    fn ids_are_stable_across_removal() {
        let mut s = NBodySystem::new(1.0);
        let a = s.add(Body::uniform("A", Material::Rock, 1.0, Vec3::ZERO, Vec3::ZERO).unwrap());
        let b = s.add(Body::uniform("B", Material::Rock, 1.0, Vec3::X, Vec3::ZERO).unwrap());
        let c = s.add(Body::uniform("C", Material::Rock, 1.0, Vec3::Y, Vec3::ZERO).unwrap());
        s.remove(b);
        assert_eq!(s.get(a).unwrap().name, "A");
        assert_eq!(s.get(c).unwrap().name, "C");
        assert!(s.get(b).is_none());
        let g0 = s.generation();
        s.find_mut("c").unwrap().name = "D".into();
        assert!(s.generation() > g0);
    }
}
