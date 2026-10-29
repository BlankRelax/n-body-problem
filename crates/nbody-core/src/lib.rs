//! # nbody-core
//!
//! Newtonian gravitational N-body dynamics with finite-size, layered bodies.
//!
//! The crate is organised bottom-up:
//!
//! | module | what it does |
//! |--------|--------------|
//! | [`vector`] | a small `f64` 3-vector |
//! | [`units`] | physical constants and unit systems (internally everything is SI) |
//! | [`material`] | bulk materials: density and coefficient of restitution |
//! | [`structure`] | concentric-shell spheres: mass, radius, moment of inertia, differentiated mergers |
//! | [`body`] | a body: name, kinematic state, structure |
//! | [`system`] | the set of bodies, `G`, softening, clock; orbit placement |
//! | [`force`] | accelerations by direct summation or a Barnes–Hut octree |
//! | [`integrator`] | symplectic (leapfrog, Yoshida, PEFRL) and RK4 time steppers |
//! | [`collision`] | merge / bounce / physical (Leinhardt–Stewart) collision handling |
//! | [`orbit`] | Keplerian elements, Kepler's equation, Hill/Roche/sphere-of-influence scales |
//! | [`diagnostics`] | energy, momentum and angular-momentum conservation tracking |
//! | [`projectile`] | hyperbolic injection and capture / ejection / merger classification |
//! | [`simulation`] | the driver tying integrator, forces, collisions and trackers together |
//! | [`scenario`] | TOML scenario files |
//! | [`presets`] | ready-made systems |
//!
//! ## A two-line simulation
//!
//! ```
//! use nbody_core::presets;
//! let built = presets::preset("solar-system").unwrap().build().unwrap();
//! let mut sim = built.into_simulation();
//! sim.advance(nbody_core::units::YEAR);
//! assert!(sim.energy_error().abs() < 1e-8);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod body;
pub mod collision;
pub mod diagnostics;
pub mod force;
pub mod integrator;
pub mod material;
pub mod orbit;
pub mod presets;
pub mod projectile;
pub mod scenario;
pub mod simulation;
pub mod structure;
pub mod system;
pub mod units;
pub mod vector;

pub use body::{Body, BodyId, BodyKind};
pub use collision::{CollisionEvent, CollisionMode, CollisionOutcome};
pub use diagnostics::{ConservationTracker, Diagnostics};
pub use force::ForceMethod;
pub use integrator::IntegratorKind;
pub use material::Material;
pub use orbit::{Elements, OrbitalElements};
pub use projectile::{Outcome, OutcomeClass, ProjectileSpec, ProjectileTracker};
pub use scenario::{Built, Scenario};
pub use simulation::{Simulation, SimulationConfig, TimestepPolicy};
pub use structure::{Layer, SphericalStructure};
pub use system::NBodySystem;
pub use units::UnitSystem;
pub use vector::Vec3;
