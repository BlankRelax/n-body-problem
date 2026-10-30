//! Scenario files: a TOML description of units, numerics, bodies and an
//! optional projectile, convertible to a ready-to-run [`Simulation`].
//!
//! ```toml
//! name = "Sun, Jupiter and a rogue planet"
//!
//! [units]
//! length = "au"
//! mass = "solar"
//! time = "year"
//!
//! [simulation]
//! integrator = "yoshida4"      # leapfrog | yoshida4 | pefrl | rk4
//! dt = 0.001                   # or [simulation.adaptive] eta = 0.01, dt_min = 1e-6, dt_max = 0.01
//! collisions = "physical"      # none | merge | bounce | physical
//! force = "direct"             # direct | barnes_hut (with theta = 0.5)
//!
//! [[bodies]]
//! name = "Sun"
//! kind = "star"
//! mass = 1.0
//! material = "plasma"
//!
//! [[bodies]]
//! name = "Jupiter"
//! mass = 9.546e-4
//! radius = 4.673e-4
//! material = "gas"             # mantle material; with `core = "rock"` a two-layer fit is made
//! core = "rock"
//! orbit = { around = "Sun", a = 5.2, e = 0.0489, inclination = 1.3, mean_anomaly = 20.0 }
//!
//! [projectile]
//! name = "Rogue"
//! mass = 3.0e-6
//! material = "rock"
//! v_infinity = 1.0             # length/time units
//! impact_parameter = 5.0       # length units
//! direction = [1, 0, 0]
//! max_time = 50.0
//! ```
//!
//! Angles are in degrees, densities of custom materials always in kg/m³,
//! everything else in the declared units. Bodies may be positioned either
//! explicitly (`position`, `velocity`) or by Keplerian elements relative to
//! a previously listed body (`orbit`).

use crate::body::{Body, BodyKind};
use crate::collision::CollisionMode;
use crate::force::{self, ForceMethod};
use crate::integrator::IntegratorKind;
use crate::material::Material;
use crate::orbit::{self, Elements};
use crate::projectile::ProjectileSpec;
use crate::simulation::{Simulation, SimulationConfig, TimestepPolicy};
use crate::structure::{self, Layer, SphericalStructure};
use crate::system::NBodySystem;
use crate::units::{UnitError, UnitSystem, G_SI};
use crate::vector::Vec3;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ScenarioError {
    #[error("TOML parse error: {0}")]
    Parse(String),
    #[error("TOML write error: {0}")]
    Serialize(String),
    #[error("I/O error: {0}")]
    Io(String),
    #[error(transparent)]
    Units(#[from] UnitError),
    #[error("body '{name}': {reason}")]
    Body { name: String, reason: String },
    #[error("unknown force method '{0}' (direct | barnes_hut)")]
    Force(String),
    #[error("projectile: {0}")]
    Projectile(String),
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct UnitSpec {
    pub length: String,
    pub mass: String,
    pub time: String,
}

impl Default for UnitSpec {
    fn default() -> Self {
        UnitSpec { length: "m".into(), mass: "kg".into(), time: "s".into() }
    }
}

impl UnitSpec {
    pub fn new(length: &str, mass: &str, time: &str) -> Self {
        UnitSpec { length: length.into(), mass: mass.into(), time: time.into() }
    }

    pub fn astronomical() -> Self {
        Self::new("au", "solar", "year")
    }

    pub fn resolve(&self) -> Result<UnitSystem, UnitError> {
        UnitSystem::from_names(&self.length, &self.mass, &self.time)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct AdaptiveSpec {
    pub eta: f64,
    pub dt_min: f64,
    pub dt_max: f64,
}

impl Default for AdaptiveSpec {
    fn default() -> Self {
        AdaptiveSpec { eta: 0.01, dt_min: 0.0, dt_max: f64::INFINITY }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct SimulationSpec {
    pub integrator: IntegratorKind,
    /// Fixed step in time units. Omitted: 1/1000 of the shortest orbital period.
    pub dt: Option<f64>,
    pub adaptive: Option<AdaptiveSpec>,
    /// Plummer softening in length units.
    pub softening: f64,
    pub collisions: CollisionMode,
    /// `direct` or `barnes_hut`.
    pub force: String,
    pub theta: f64,
    /// Override of G in the scenario's units (e.g. `1.0` for textbook problems).
    pub g: Option<f64>,
    /// Suggested run length in time units.
    pub duration: Option<f64>,
    /// Shift to the barycentric frame after building.
    pub barycentric: bool,
}

impl Default for SimulationSpec {
    fn default() -> Self {
        SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: None,
            adaptive: None,
            softening: 0.0,
            collisions: CollisionMode::Physical,
            force: "direct".into(),
            theta: 0.5,
            g: None,
            duration: None,
            barycentric: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LayerSpec {
    pub outer_radius: f64,
    pub material: Material,
}

/// Keplerian placement relative to another body; angles in degrees.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct OrbitSpec {
    pub around: String,
    pub a: f64,
    pub e: f64,
    pub inclination: f64,
    pub raan: f64,
    pub arg_periapsis: f64,
    pub true_anomaly: Option<f64>,
    pub mean_anomaly: Option<f64>,
}

impl Default for OrbitSpec {
    fn default() -> Self {
        OrbitSpec { around: String::new(), a: 1.0, e: 0.0, inclination: 0.0, raan: 0.0, arg_periapsis: 0.0, true_anomaly: None, mean_anomaly: None }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct BodySpec {
    pub name: String,
    pub kind: Option<BodyKind>,
    pub mass: Option<f64>,
    pub radius: Option<f64>,
    /// Material of a uniform body, or of the mantle when `core` is set.
    pub material: Option<Material>,
    /// With `mass`, `radius` and `material`: fit a core of this material.
    pub core: Option<Material>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub layers: Vec<LayerSpec>,
    pub position: Option<Vec3>,
    pub velocity: Option<Vec3>,
    pub orbit: Option<OrbitSpec>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ProjectileFileSpec {
    pub name: String,
    pub mass: Option<f64>,
    pub radius: Option<f64>,
    pub material: Option<Material>,
    pub core: Option<Material>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub layers: Vec<LayerSpec>,
    pub v_infinity: f64,
    pub impact_parameter: f64,
    pub direction: Vec3,
    pub offset_direction: Option<Vec3>,
    pub start_distance: Option<f64>,
    /// How long to run before giving up on a verdict, time units.
    pub max_time: Option<f64>,
}

impl Default for ProjectileFileSpec {
    fn default() -> Self {
        ProjectileFileSpec {
            name: "projectile".into(),
            mass: None,
            radius: None,
            material: None,
            core: None,
            layers: Vec::new(),
            v_infinity: 1.0,
            impact_parameter: 0.0,
            direction: Vec3::X,
            offset_direction: None,
            start_distance: None,
            max_time: None,
        }
    }
}

/// Hints for an interactive viewer. The physics ignores them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ViewSpec {
    /// Draw the Lagrange points of this `[primary, secondary]` pair of bodies.
    pub lagrange_points: Option<[String; 2]>,
    /// Start in the frame that co-rotates with that pair.
    pub corotating: bool,
    /// Simulated time the trails should cover, in time units (default: the run's duration).
    pub trail_span: Option<f64>,
    /// Initial tilt of the view out of the x-y plane, degrees.
    pub tilt: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct Scenario {
    pub name: String,
    pub description: String,
    pub units: UnitSpec,
    pub simulation: SimulationSpec,
    pub bodies: Vec<BodySpec>,
    pub projectile: Option<ProjectileFileSpec>,
    pub view: ViewSpec,
}

/// Everything a scenario resolves to, in SI.
#[derive(Debug)]
pub struct Built {
    pub system: NBodySystem,
    pub config: SimulationConfig,
    pub units: UnitSystem,
    pub unit_spec: UnitSpec,
    pub projectile: Option<ProjectileSpec>,
    pub projectile_max_time: Option<f64>,
    pub duration: Option<f64>,
}

impl Built {
    pub fn into_simulation(self) -> Simulation {
        Simulation::new(self.system, self.config)
    }
}

fn build_structure(
    label: &str,
    mass: Option<f64>,
    radius: Option<f64>,
    material: &Option<Material>,
    core: &Option<Material>,
    layers: &[LayerSpec],
    u: &UnitSystem,
) -> Result<SphericalStructure, ScenarioError> {
    let err = |reason: String| ScenarioError::Body { name: label.to_string(), reason };
    if !layers.is_empty() {
        let ls = layers.iter().map(|l| Layer { outer_radius: l.outer_radius * u.length, material: l.material.clone() }).collect();
        return SphericalStructure::from_layers(ls).map_err(|e| err(e.to_string()));
    }
    let mass_si = mass.map(|m| m * u.mass);
    let radius_si = radius.map(|r| r * u.length);
    match (mass_si, radius_si) {
        (Some(m), None) => SphericalStructure::uniform_with_mass(material.clone().unwrap_or(Material::Rock), m).map_err(|e| err(e.to_string())),
        (None, Some(r)) => SphericalStructure::uniform_with_radius(material.clone().unwrap_or(Material::Rock), r).map_err(|e| err(e.to_string())),
        (Some(m), Some(r)) => {
            if let (Some(core), Some(mantle)) = (core, material) {
                return structure::fit_two_layer(core.clone(), mantle.clone(), m, r).map_err(|e| err(e.to_string()));
            }
            let density = m / structure::sphere_volume(r);
            let base = material.clone().unwrap_or(Material::Rock);
            let mut custom = Material::custom(base.name(), density, base.restitution());
            if let Material::Custom { colour, .. } = &mut custom {
                *colour = base.colour();
            }
            SphericalStructure::uniform_with_radius(custom, r).map_err(|e| err(e.to_string()))
        }
        (None, None) => Err(err("needs a mass, a radius or layers".into())),
    }
}

fn parse_force(name: &str, theta: f64) -> Result<ForceMethod, ScenarioError> {
    match name.trim().to_ascii_lowercase().as_str() {
        "direct" | "n2" | "exact" => Ok(ForceMethod::Direct),
        "barnes_hut" | "barnes-hut" | "bh" | "tree" => Ok(ForceMethod::BarnesHut { theta }),
        other => Err(ScenarioError::Force(other.to_string())),
    }
}

impl Scenario {
    pub fn from_toml(text: &str) -> Result<Scenario, ScenarioError> {
        toml::from_str(text).map_err(|e| ScenarioError::Parse(e.to_string()))
    }

    pub fn to_toml(&self) -> Result<String, ScenarioError> {
        toml::to_string_pretty(self).map_err(|e| ScenarioError::Serialize(e.to_string()))
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Scenario, ScenarioError> {
        let text = std::fs::read_to_string(path.as_ref()).map_err(|e| ScenarioError::Io(format!("{}: {e}", path.as_ref().display())))?;
        Self::from_toml(&text)
    }

    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), ScenarioError> {
        std::fs::write(path.as_ref(), self.to_toml()?).map_err(|e| ScenarioError::Io(format!("{}: {e}", path.as_ref().display())))
    }

    /// Resolve units, build the system and the numerical configuration.
    pub fn build(&self) -> Result<Built, ScenarioError> {
        let u = self.units.resolve()?;
        let g = match self.simulation.g {
            Some(gu) => gu * u.length.powi(3) / (u.mass * u.time * u.time),
            None => G_SI,
        };
        let mut system = NBodySystem::new(g);
        system.softening = self.simulation.softening * u.length;

        for spec in &self.bodies {
            let name = if spec.name.is_empty() { format!("body{}", system.len() + 1) } else { spec.name.clone() };
            let err = |reason: String| ScenarioError::Body { name: name.clone(), reason };
            if system.find(&name).is_some() {
                return Err(err("duplicate name".into()));
            }
            let structure = build_structure(&name, spec.mass, spec.radius, &spec.material, &spec.core, &spec.layers, &u)?;
            let (position, velocity) = if let Some(o) = &spec.orbit {
                let primary = system.find(&o.around).ok_or_else(|| err(format!("orbit.around = '{}' is not a body listed above it", o.around)))?.clone();
                let mu = g * (primary.mass() + structure.mass());
                let nu = match (o.true_anomaly, o.mean_anomaly) {
                    (Some(t), _) => t.to_radians(),
                    (None, Some(m)) => orbit::mean_to_true_anomaly(m.to_radians(), o.e),
                    (None, None) => 0.0,
                };
                let el = Elements {
                    a: o.a * u.length,
                    e: o.e,
                    inclination: o.inclination.to_radians(),
                    raan: o.raan.to_radians(),
                    arg_periapsis: o.arg_periapsis.to_radians(),
                    true_anomaly: nu,
                };
                let (r, v) = orbit::state_from_elements(mu, &el).map_err(|e| err(e.to_string()))?;
                (primary.position + r, primary.velocity + v)
            } else {
                (spec.position.unwrap_or(Vec3::ZERO) * u.length, spec.velocity.unwrap_or(Vec3::ZERO) * u.velocity())
            };
            let body = Body::new(name.clone(), structure, position, velocity).with_kind(spec.kind.unwrap_or_default());
            system.add(body);
        }
        if self.simulation.barycentric && !system.is_empty() {
            system.move_to_barycentre();
        }

        let timestep = if let Some(a) = self.simulation.adaptive {
            TimestepPolicy::Adaptive { eta: a.eta, dt_min: a.dt_min * u.time, dt_max: a.dt_max * u.time }
        } else {
            let dt = match self.simulation.dt {
                Some(dt) => dt * u.time,
                None => {
                    let t = force::min_dynamical_time(&system.positions(), &system.masses(), g, system.softening);
                    if t.is_finite() {
                        2.0 * std::f64::consts::PI * t / 1000.0
                    } else {
                        u.time
                    }
                }
            };
            TimestepPolicy::Fixed { dt }
        };
        let config = SimulationConfig {
            integrator: self.simulation.integrator,
            timestep,
            force: parse_force(&self.simulation.force, self.simulation.theta)?,
            collisions: self.simulation.collisions,
        };

        let (projectile, projectile_max_time) = match &self.projectile {
            Some(p) => {
                let structure = build_structure(&p.name, p.mass, p.radius, &p.material, &p.core, &p.layers, &u)
                    .map_err(|e| ScenarioError::Projectile(e.to_string()))?;
                let spec = ProjectileSpec {
                    name: p.name.clone(),
                    structure,
                    v_infinity: p.v_infinity * u.velocity(),
                    impact_parameter: p.impact_parameter * u.length,
                    direction: p.direction,
                    offset_direction: p.offset_direction,
                    start_distance: p.start_distance.map(|d| d * u.length),
                };
                (Some(spec), p.max_time.map(|t| t * u.time))
            }
            None => (None, None),
        };

        Ok(Built {
            system,
            config,
            units: u,
            unit_spec: self.units.clone(),
            projectile,
            projectile_max_time,
            duration: self.simulation.duration.map(|d| d * u.time),
        })
    }

    /// Snapshot a live system as a scenario with explicit states and layers.
    pub fn snapshot(system: &NBodySystem, config: &SimulationConfig, unit_spec: &UnitSpec, name: &str) -> Result<Scenario, ScenarioError> {
        let u = unit_spec.resolve()?;
        let bodies = system
            .bodies()
            .iter()
            .map(|b| BodySpec {
                name: b.name.clone(),
                kind: Some(b.kind),
                layers: b
                    .structure()
                    .layers()
                    .iter()
                    .map(|l| LayerSpec { outer_radius: l.outer_radius / u.length, material: l.material.clone() })
                    .collect(),
                position: Some(b.position / u.length),
                velocity: Some(b.velocity / u.velocity()),
                ..Default::default()
            })
            .collect();
        let (dt, adaptive) = match config.timestep {
            TimestepPolicy::Fixed { dt } => (Some(dt / u.time), None),
            TimestepPolicy::Adaptive { eta, dt_min, dt_max } => (None, Some(AdaptiveSpec { eta, dt_min: dt_min / u.time, dt_max: dt_max / u.time })),
        };
        let (force, theta) = match config.force {
            ForceMethod::Direct => ("direct".to_string(), 0.5),
            ForceMethod::BarnesHut { theta } => ("barnes_hut".to_string(), theta),
        };
        let g_units = system.g * u.mass * u.time * u.time / u.length.powi(3);
        let g = if (system.g / G_SI - 1.0).abs() > 1e-9 { Some(g_units) } else { None };
        Ok(Scenario {
            name: name.to_string(),
            description: format!("snapshot at t = {} s", system.time),
            units: unit_spec.clone(),
            simulation: SimulationSpec {
                integrator: config.integrator,
                dt,
                adaptive,
                softening: system.softening / u.length,
                collisions: config.collisions,
                force,
                theta,
                g,
                duration: None,
                barycentric: false,
            },
            bodies,
            projectile: None,
            view: ViewSpec::default(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::*;
    use approx::assert_relative_eq;

    const SAMPLE: &str = r#"
name = "test"
[units]
length = "au"
mass = "solar"
time = "year"

[simulation]
integrator = "pefrl"
dt = 0.001
collisions = "merge"

[[bodies]]
name = "Sun"
kind = "star"
mass = 1.0
material = "plasma"

[[bodies]]
name = "Earth"
mass = 3.003e-6
radius = 4.2635e-5
material = "rock"
core = "iron"
orbit = { around = "Sun", a = 1.0, e = 0.0167, mean_anomaly = 90.0 }

[projectile]
name = "Rogue"
mass = 1.0e-6
material = "ice"
v_infinity = 2.0
impact_parameter = 0.5
max_time = 20
"#;

    #[test]
    fn parse_and_build_sample() {
        let sc = Scenario::from_toml(SAMPLE).unwrap();
        let built = sc.build().unwrap();
        assert_eq!(built.system.len(), 2);
        let earth = built.system.find("Earth").unwrap();
        assert_relative_eq!(earth.mass(), 3.003e-6 * M_SUN, max_relative = 1e-9);
        assert_eq!(earth.structure().layers().len(), 2, "core/mantle fit");
        assert_relative_eq!(earth.radius(), 4.2635e-5 * AU, max_relative = 1e-9);
        let sun = built.system.find("Sun").unwrap();
        let rel = earth.position - sun.position;
        let a = AU;
        assert!(rel.norm() > a * (1.0 - 0.0167) && rel.norm() < a * (1.0 + 0.0167));
        assert_relative_eq!(built.config.timestep.nominal_dt(), 0.001 * YEAR, max_relative = 1e-12);
        assert_eq!(built.config.integrator, IntegratorKind::Pefrl);
        assert_eq!(built.config.collisions, CollisionMode::Merge);
        let p = built.projectile.unwrap();
        assert_relative_eq!(p.v_infinity, 2.0 * AU / YEAR, max_relative = 1e-12);
        assert_relative_eq!(built.projectile_max_time.unwrap(), 20.0 * YEAR, max_relative = 1e-12);
        // Barycentric frame.
        assert!(built.system.momentum().norm() < 1e-3 * M_SUN);
    }

    #[test]
    fn snapshot_round_trip_preserves_state() {
        let sc = Scenario::from_toml(SAMPLE).unwrap();
        let built = sc.build().unwrap();
        let snap = Scenario::snapshot(&built.system, &built.config, &built.unit_spec, "snap").unwrap();
        let text = snap.to_toml().unwrap();
        let again = Scenario::from_toml(&text).unwrap().build().unwrap();
        for (a, b) in built.system.bodies().iter().zip(again.system.bodies()) {
            assert_eq!(a.name, b.name);
            assert!((a.position - b.position).norm() < 1e-6 * AU);
            assert!((a.velocity - b.velocity).norm() < 1e-9 * AU / YEAR);
            assert_relative_eq!(a.mass(), b.mass(), max_relative = 1e-9);
        }
    }

    #[test]
    fn g_override_and_auto_dt() {
        let text = r#"
[simulation]
g = 1.0
[[bodies]]
name = "a"
mass = 1.0
position = [0, 0, 0]
[[bodies]]
name = "b"
mass = 1.0
position = [1, 0, 0]
velocity = [0, 1, 0]
"#;
        let built = Scenario::from_toml(text).unwrap().build().unwrap();
        assert_eq!(built.system.g, 1.0);
        // t_dyn = sqrt(1 / 2) → dt = 2π t_dyn / 1000
        let expected = 2.0 * std::f64::consts::PI * (0.5f64).sqrt() / 1000.0;
        assert_relative_eq!(built.config.timestep.nominal_dt(), expected, max_relative = 1e-12);
    }

    #[test]
    fn errors_are_descriptive() {
        let text = r#"
[[bodies]]
name = "x"
orbit = { around = "nobody", a = 1 }
mass = 1
"#;
        let e = Scenario::from_toml(text).unwrap().build().unwrap_err();
        assert!(e.to_string().contains("nobody"), "{e}");
        let e = Scenario::from_toml("[units]\nlength = \"cubit\"").unwrap().build().unwrap_err();
        assert!(e.to_string().contains("cubit"));
    }
}
