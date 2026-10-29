//! Ready-made scenarios: the Solar System, textbook three-body problems,
//! a star cluster and projectile experiments.

use crate::body::BodyKind;
use crate::collision::CollisionMode;
use crate::integrator::IntegratorKind;
use crate::material::Material;
use crate::orbit;
use crate::projectile;
use crate::scenario::{AdaptiveSpec, BodySpec, OrbitSpec, ProjectileFileSpec, Scenario, SimulationSpec, UnitSpec};
use crate::units::*;
use crate::vector::Vec3;
use rand::{Rng, SeedableRng};

pub const NAMES: [(&str, &str); 9] = [
    ("solar-system", "Sun and the eight planets at J2000 (au, solar masses, years)"),
    ("inner-solar-system", "Sun, Mercury to Mars, and the Moon around Earth"),
    ("earth-moon", "Earth and Moon (km, kg, hours)"),
    ("figure-eight", "Chenciner-Montgomery figure-eight choreography, G = 1"),
    ("pythagorean", "Burrau's Pythagorean three-body problem (3-4-5 triangle), G = 1"),
    ("binary-planet", "Alpha-Centauri-like binary with an Earth-mass planet around the primary"),
    ("plummer", "200-body Plummer-sphere star cluster, Barnes-Hut forces (pc, solar masses, Myr)"),
    ("projectile-capture", "Sun + Jupiter visited by a rogue Earth-mass planet aimed at Jupiter's orbit"),
    ("projectile-flyby", "Sun + Earth visited by a fast comet"),
];

pub fn names() -> &'static [(&'static str, &'static str)] {
    &NAMES
}

pub fn preset(name: &str) -> Option<Scenario> {
    Some(match name.trim().to_ascii_lowercase().as_str() {
        "solar-system" | "solar_system" | "solar" => solar_system(),
        "inner-solar-system" | "inner" => inner_solar_system(),
        "earth-moon" | "earth_moon" => earth_moon(),
        "figure-eight" | "figure8" | "figure_eight" => figure_eight(),
        "pythagorean" | "burrau" => pythagorean(),
        "binary-planet" | "binary" => binary_planet(),
        "plummer" | "cluster" => plummer(200, 7),
        "projectile-capture" | "capture" => projectile_capture(),
        "projectile-flyby" | "flyby" => projectile_flyby(),
        _ => return None,
    })
}

struct Planet {
    name: &'static str,
    mass: f64,
    radius: f64,
    mantle: Material,
    core: Material,
    // JPL approximate Keplerian elements at J2000 (au, degrees).
    a: f64,
    e: f64,
    i: f64,
    long_node: f64,
    long_peri: f64,
    mean_long: f64,
}

const PLANETS: [Planet; 8] = [
    Planet { name: "Mercury", mass: 3.3011e23, radius: 2.4397e6, mantle: Material::Rock, core: Material::Iron, a: 0.38709927, e: 0.20563593, i: 7.00497902, long_node: 48.33076593, long_peri: 77.45779628, mean_long: 252.25032350 },
    Planet { name: "Venus", mass: 4.8675e24, radius: 6.0518e6, mantle: Material::Rock, core: Material::Iron, a: 0.72333566, e: 0.00677672, i: 3.39467605, long_node: 76.67984255, long_peri: 131.60246718, mean_long: 181.97909950 },
    Planet { name: "Earth", mass: 5.9724e24, radius: 6.3710e6, mantle: Material::Rock, core: Material::Iron, a: 1.00000261, e: 0.01671123, i: -0.00001531, long_node: 0.0, long_peri: 102.93768193, mean_long: 100.46457166 },
    Planet { name: "Mars", mass: 6.4171e23, radius: 3.3895e6, mantle: Material::Rock, core: Material::Iron, a: 1.52371034, e: 0.09339410, i: 1.84969142, long_node: 49.55953891, long_peri: -23.94362959, mean_long: -4.55343205 },
    Planet { name: "Jupiter", mass: 1.8982e27, radius: 6.9911e7, mantle: Material::Gas, core: Material::Rock, a: 5.20288700, e: 0.04838624, i: 1.30439695, long_node: 100.47390909, long_peri: 14.72847983, mean_long: 34.39644051 },
    Planet { name: "Saturn", mass: 5.6834e26, radius: 5.8232e7, mantle: Material::Gas, core: Material::Rock, a: 9.53667594, e: 0.05386179, i: 2.48599187, long_node: 113.66242448, long_peri: 92.59887831, mean_long: 49.95424423 },
    Planet { name: "Uranus", mass: 8.6810e25, radius: 2.5362e7, mantle: Material::Gas, core: Material::Ice, a: 19.18916464, e: 0.04725744, i: 0.77263783, long_node: 74.01692503, long_peri: 170.95427630, mean_long: 313.23810451 },
    Planet { name: "Neptune", mass: 1.02413e26, radius: 2.4622e7, mantle: Material::Gas, core: Material::Ice, a: 30.06992276, e: 0.00859048, i: 1.77004347, long_node: 131.78422574, long_peri: 44.96476227, mean_long: -55.12002969 },
];

fn planet_spec(p: &Planet) -> BodySpec {
    BodySpec {
        name: p.name.into(),
        kind: Some(BodyKind::Planet),
        mass: Some(p.mass / M_SUN),
        radius: Some(p.radius / AU),
        material: Some(p.mantle.clone()),
        core: Some(p.core.clone()),
        orbit: Some(OrbitSpec {
            around: "Sun".into(),
            a: p.a,
            e: p.e,
            inclination: p.i,
            raan: p.long_node,
            arg_periapsis: p.long_peri - p.long_node,
            mean_anomaly: Some(p.mean_long - p.long_peri),
            true_anomaly: None,
        }),
        ..Default::default()
    }
}

fn sun_spec() -> BodySpec {
    BodySpec {
        name: "Sun".into(),
        kind: Some(BodyKind::Star),
        mass: Some(1.0),
        radius: Some(R_SUN / AU),
        material: Some(Material::Plasma),
        ..Default::default()
    }
}

fn moon_spec() -> BodySpec {
    BodySpec {
        name: "Moon".into(),
        kind: Some(BodyKind::Moon),
        mass: Some(M_MOON / M_SUN),
        radius: Some(R_MOON / AU),
        material: Some(Material::Rock),
        core: Some(Material::Iron),
        orbit: Some(OrbitSpec {
            around: "Earth".into(),
            a: 384_400e3 / AU,
            e: 0.0549,
            inclination: 5.145,
            raan: 125.08,
            arg_periapsis: 318.15,
            mean_anomaly: Some(135.27),
            true_anomaly: None,
        }),
        ..Default::default()
    }
}

pub fn solar_system() -> Scenario {
    Scenario {
        name: "Solar System".into(),
        description: "Sun and the eight planets from JPL's J2000 mean elements. Two-layer planets (core + mantle) reproduce each body's bulk mass and radius.".into(),
        units: UnitSpec::astronomical(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: Some(0.0005),
            duration: Some(12.0),
            ..Default::default()
        },
        bodies: std::iter::once(sun_spec()).chain(PLANETS.iter().map(planet_spec)).collect(),
        projectile: None,
    }
}

pub fn inner_solar_system() -> Scenario {
    let mut bodies = vec![sun_spec()];
    bodies.extend(PLANETS.iter().take(4).map(planet_spec));
    bodies.push(moon_spec());
    Scenario {
        name: "Inner Solar System".into(),
        description: "Sun, Mercury, Venus, Earth (with the Moon) and Mars.".into(),
        units: UnitSpec::astronomical(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: Some(0.0001),
            duration: Some(2.0),
            ..Default::default()
        },
        bodies,
        projectile: None,
    }
}

pub fn earth_moon() -> Scenario {
    Scenario {
        name: "Earth and Moon".into(),
        description: "The Earth-Moon system in km, kg and hours; one sidereal month is 655.7 h.".into(),
        units: UnitSpec::new("km", "kg", "hour"),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: Some(0.05),
            duration: Some(655.7),
            ..Default::default()
        },
        bodies: vec![
            BodySpec {
                name: "Earth".into(),
                kind: Some(BodyKind::Planet),
                mass: Some(M_EARTH),
                radius: Some(R_EARTH / 1e3),
                material: Some(Material::Rock),
                core: Some(Material::Iron),
                ..Default::default()
            },
            BodySpec {
                name: "Moon".into(),
                kind: Some(BodyKind::Moon),
                mass: Some(M_MOON),
                radius: Some(R_MOON / 1e3),
                material: Some(Material::Rock),
                core: Some(Material::Iron),
                orbit: Some(OrbitSpec { around: "Earth".into(), a: 384_400.0, e: 0.0549, inclination: 5.145, ..Default::default() }),
                ..Default::default()
            },
        ],
        projectile: None,
    }
}

fn unit_body(name: &str, mass: f64, position: [f64; 3], velocity: [f64; 3]) -> BodySpec {
    BodySpec {
        name: name.into(),
        kind: Some(BodyKind::Star),
        mass: Some(mass),
        material: Some(Material::Plasma),
        position: Some(Vec3::from(position)),
        velocity: Some(Vec3::from(velocity)),
        ..Default::default()
    }
}

/// Chenciner & Montgomery (2000); initial conditions from Šuvakov & Dmitrašinović.
pub fn figure_eight() -> Scenario {
    let x = 0.97000436;
    let y = 0.24308753;
    let vx = 0.93240737;
    let vy = 0.86473146;
    Scenario {
        name: "Figure eight".into(),
        description: "Three equal masses chasing each other around a figure-eight (period 6.3259), G = 1.".into(),
        units: UnitSpec::default(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            dt: Some(0.002),
            g: Some(1.0),
            duration: Some(3.0 * 6.32591398),
            collisions: CollisionMode::None,
            ..Default::default()
        },
        bodies: vec![
            unit_body("A", 1.0, [-x, y, 0.0], [vx / 2.0, vy / 2.0, 0.0]),
            unit_body("B", 1.0, [x, -y, 0.0], [vx / 2.0, vy / 2.0, 0.0]),
            unit_body("C", 1.0, [0.0, 0.0, 0.0], [-vx, -vy, 0.0]),
        ],
        projectile: None,
    }
}

/// Burrau (1913); Szebehely & Peters (1967). Ends with a binary and an escaper near t ≈ 60.
pub fn pythagorean() -> Scenario {
    Scenario {
        name: "Pythagorean problem".into(),
        description: "Masses 3, 4, 5 at rest at the vertices of a 3-4-5 right triangle, G = 1. Chaotic close encounters until, near t = 60, the lightest body (mass 3) is ejected and masses 4 and 5 recoil as a tight binary.".into(),
        units: UnitSpec::default(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            adaptive: Some(AdaptiveSpec { eta: 0.003, dt_min: 1e-10, dt_max: 0.01 }),
            g: Some(1.0),
            duration: Some(70.0),
            collisions: CollisionMode::None,
            ..Default::default()
        },
        bodies: vec![
            unit_body("m3", 3.0, [1.0, 3.0, 0.0], [0.0; 3]),
            unit_body("m4", 4.0, [-2.0, -1.0, 0.0], [0.0; 3]),
            unit_body("m5", 5.0, [1.0, -1.0, 0.0], [0.0; 3]),
        ],
        projectile: None,
    }
}

pub fn binary_planet() -> Scenario {
    Scenario {
        name: "Binary star with a planet".into(),
        description: "Alpha-Centauri-like pair (1.1 + 0.9 solar masses, a = 23.4 au, e = 0.52) with an Earth-mass planet 1 au from the primary; S-type orbits here are stable for billions of years.".into(),
        units: UnitSpec::astronomical(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: Some(0.001),
            duration: Some(160.0),
            ..Default::default()
        },
        bodies: vec![
            BodySpec {
                name: "Alpha A".into(),
                kind: Some(BodyKind::Star),
                mass: Some(1.1),
                radius: Some(1.22 * R_SUN / AU),
                material: Some(Material::Plasma),
                ..Default::default()
            },
            BodySpec {
                name: "Alpha B".into(),
                kind: Some(BodyKind::Star),
                mass: Some(0.9),
                radius: Some(0.86 * R_SUN / AU),
                material: Some(Material::Plasma),
                orbit: Some(OrbitSpec { around: "Alpha A".into(), a: 23.4, e: 0.52, inclination: 0.0, mean_anomaly: Some(180.0), ..Default::default() }),
                ..Default::default()
            },
            BodySpec {
                name: "Planet".into(),
                kind: Some(BodyKind::Planet),
                mass: Some(M_EARTH / M_SUN),
                radius: Some(R_EARTH / AU),
                material: Some(Material::Rock),
                core: Some(Material::Iron),
                orbit: Some(OrbitSpec { around: "Alpha A".into(), a: 1.0, e: 0.02, ..Default::default() }),
                ..Default::default()
            },
        ],
        projectile: None,
    }
}

/// Plummer (1911) sphere sampled with the Aarseth, Hénon & Wielen (1974) recipe.
pub fn plummer(n: usize, seed: u64) -> Scenario {
    let units = UnitSpec::new("pc", "solar", "myr");
    let u = units.resolve().expect("known units");
    let total_mass = 5.0 * n as f64; // 5 solar masses each
    let scale_radius = 1.0; // pc
    let g = u.g();
    let v_scale = (g * total_mass / scale_radius).sqrt();
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut bodies = Vec::with_capacity(n);
    let random_direction = |rng: &mut rand::rngs::StdRng, r: f64| {
        let z = (1.0 - 2.0 * rng.gen::<f64>()) * r;
        let phi = 2.0 * std::f64::consts::PI * rng.gen::<f64>();
        let rho = (r * r - z * z).max(0.0).sqrt();
        Vec3::new(rho * phi.cos(), rho * phi.sin(), z)
    };
    for k in 0..n {
        let r = loop {
            let x1: f64 = rng.gen_range(1e-6..1.0);
            let r = (x1.powf(-2.0 / 3.0) - 1.0).powf(-0.5);
            if r < 10.0 {
                break r;
            }
        };
        let q = loop {
            let x4: f64 = rng.gen();
            let x5: f64 = rng.gen();
            if 0.1 * x5 < x4 * x4 * (1.0 - x4 * x4).powf(3.5) {
                break x4;
            }
        };
        let v = q * 2f64.sqrt() * (1.0 + r * r).powf(-0.25);
        let position = random_direction(&mut rng, r) * scale_radius;
        let velocity = random_direction(&mut rng, v) * v_scale;
        bodies.push(BodySpec {
            name: format!("star{:03}", k + 1),
            kind: Some(BodyKind::Star),
            mass: Some(total_mass / n as f64),
            material: Some(Material::Plasma),
            position: Some(position),
            velocity: Some(velocity),
            ..Default::default()
        });
    }
    Scenario {
        name: format!("Plummer sphere, N = {n}"),
        description: "A self-consistent star cluster in virial equilibrium (2K/|U| ≈ 1). Uses the Barnes-Hut tree with softening to tame close pairs.".into(),
        units,
        simulation: SimulationSpec {
            integrator: IntegratorKind::Leapfrog,
            dt: Some(0.002),
            softening: 0.02,
            force: "barnes_hut".into(),
            theta: 0.5,
            collisions: CollisionMode::None,
            duration: Some(10.0),
            ..Default::default()
        },
        bodies,
        projectile: None,
    }
}

/// Jupiter's phase offset (degrees of orbital longitude) from the point
/// where the rogue planet reaches periapsis, at the moment it gets there.
/// Chosen by scanning (`scan_capture_lag`) so the default run is a capture.
pub const CAPTURE_LAG_DEG: f64 = -2.5;

pub fn projectile_capture() -> Scenario {
    projectile_capture_with_lag(CAPTURE_LAG_DEG)
}

/// The rogue-planet scenario with Jupiter placed `lag_deg` ahead of the
/// projectile's periapsis point at the analytically predicted arrival time.
///
/// The approach hyperbola is aimed so that its periapsis sits exactly on
/// Jupiter's orbit (`b = q sqrt(1 + 2μ/(q v∞²))`), and the time of flight to
/// periapsis follows from the hyperbolic Kepler equation, so the encounter
/// geometry is fixed by two numbers: the lag angle and the excess speed.
pub fn projectile_capture_with_lag(lag_deg: f64) -> Scenario {
    let jupiter = &PLANETS[4];
    let units = UnitSpec::astronomical();
    let u = units.resolve().expect("known units");
    let g = u.g();
    let (m_sun, m_j, m_p) = (1.0, jupiter.mass / M_SUN, M_EARTH / M_SUN);
    let a_j = 5.2029;
    let (v_inf, start) = (1.0, 60.0);
    let mu = g * (m_sun + m_j + m_p);
    let b = orbit::impact_parameter_for_periapsis(mu, v_inf, a_j);
    let (dir, off) = (Vec3::X, Vec3::new(0.0, -1.0, 0.0));
    let geo = projectile::approach_geometry(mu, v_inf, b, start, dir, Some(off)).expect("valid approach");
    let phi = geo.periapsis.y.atan2(geo.periapsis.x);
    let n_j = (g * (m_sun + m_j) / (a_j * a_j * a_j)).sqrt();
    let nu_j = (phi - n_j * geo.time_to_periapsis + lag_deg.to_radians()).rem_euclid(2.0 * std::f64::consts::PI).to_degrees();
    let mut jup = planet_spec(jupiter);
    jup.orbit = Some(OrbitSpec { around: "Sun".into(), a: a_j, e: 0.0, true_anomaly: Some(nu_j), ..Default::default() });
    Scenario {
        name: "Rogue planet meets Jupiter".into(),
        description: format!(
            "A rogue Earth-mass planet arrives at 1 au/yr (4.7 km/s) on a hyperbola whose periapsis lies on Jupiter's orbit; Jupiter is phased to be {lag_deg:+.1}° from that point when the rogue arrives (after {:.2} yr). Jupiter can brake it into a bound orbit, fling it out faster, or swallow it: vary impact_parameter and v_infinity (or run `nbody sweep`) to map the regimes.",
            geo.time_to_periapsis
        ),
        units,
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            adaptive: Some(AdaptiveSpec { eta: 0.01, dt_min: 1e-8, dt_max: 0.002 }),
            duration: Some(120.0),
            ..Default::default()
        },
        bodies: vec![sun_spec(), jup],
        projectile: Some(ProjectileFileSpec {
            name: "Rogue".into(),
            mass: Some(m_p),
            radius: Some(R_EARTH / AU),
            material: Some(Material::Rock),
            core: Some(Material::Iron),
            v_infinity: v_inf,
            impact_parameter: b,
            direction: dir,
            offset_direction: Some(off),
            start_distance: Some(start),
            max_time: Some(120.0),
            ..Default::default()
        }),
    }
}

pub fn projectile_flyby() -> Scenario {
    let earth = planet_spec(&PLANETS[2]);
    Scenario {
        name: "Comet flyby".into(),
        description: "A 10-km icy comet nucleus passes the Sun at 3 au/yr with a 0.3 au impact parameter: a pure hyperbolic flyby, deflected by the Rutherford angle 2 atan(GM / b v^2).".into(),
        units: UnitSpec::astronomical(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            adaptive: Some(AdaptiveSpec { eta: 0.01, dt_min: 1e-8, dt_max: 0.002 }),
            duration: Some(30.0),
            ..Default::default()
        },
        bodies: vec![sun_spec(), earth],
        projectile: Some(ProjectileFileSpec {
            name: "Comet".into(),
            radius: Some(5e3 / AU),
            material: Some(Material::Ice),
            v_infinity: 3.0,
            impact_parameter: 0.3,
            direction: Vec3::X,
            start_distance: Some(30.0),
            max_time: Some(30.0),
            ..Default::default()
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::Diagnostics;

    #[test]
    fn all_presets_build() {
        for (name, _) in NAMES {
            let sc = preset(name).unwrap_or_else(|| panic!("missing preset {name}"));
            let built = sc.build().unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(!built.system.is_empty());
            let text = sc.to_toml().unwrap();
            Scenario::from_toml(&text).unwrap().build().unwrap_or_else(|e| panic!("{name} round trip: {e}"));
        }
    }

    #[test]
    fn solar_system_periods_match_kepler() {
        let built = solar_system().build().unwrap();
        let sun = built.system.find("Sun").unwrap();
        for p in &PLANETS {
            let b = built.system.find(p.name).unwrap();
            let mu = built.system.g * (sun.mass() + b.mass());
            let el = crate::orbit::elements_from_state(mu, b.position - sun.position, b.velocity - sun.velocity);
            let period_yr = el.period.unwrap() / YEAR;
            let expected = p.a.powf(1.5);
            assert!((period_yr / expected - 1.0).abs() < 2e-3, "{}: {period_yr} vs {expected}", p.name);
            assert!((el.elements.e - p.e).abs() < 1e-6);
        }
    }

    #[test]
    fn plummer_is_near_virial_equilibrium() {
        let built = plummer(400, 3).build().unwrap();
        let d = Diagnostics::compute(&built.system);
        assert!((d.virial_ratio - 1.0).abs() < 0.25, "virial ratio {}", d.virial_ratio);
        assert!(d.total < 0.0);
    }

    fn capture_outcome(lag: f64) -> crate::projectile::OutcomeClass {
        let built = projectile_capture_with_lag(lag).build().unwrap();
        let spec = built.projectile.clone().unwrap();
        let max_t = built.projectile_max_time.unwrap();
        let mut sim = built.into_simulation();
        sim.inject_projectile(&spec).unwrap();
        sim.run_projectile(max_t).class()
    }

    /// The three regimes of the Jupiter encounter: a pass in front brakes the
    /// rogue planet into a bound orbit, a direct hit merges, a pass behind
    /// slingshots it out faster than it came.
    #[test]
    fn jupiter_encounter_regimes() {
        use crate::projectile::OutcomeClass;
        assert_eq!(capture_outcome(CAPTURE_LAG_DEG), OutcomeClass::Captured);
        assert_eq!(capture_outcome(0.0), OutcomeClass::Merged);
        assert_eq!(capture_outcome(4.5), OutcomeClass::Ejected);
    }

    /// Scan Jupiter's phase lag and print the outcome of each encounter.
    /// Run with `cargo test -p nbody-core scan_capture_lag -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn scan_capture_lag() {
        for k in -60..=60 {
            let lag = k as f64 * 0.1;
            let built = projectile_capture_with_lag(lag).build().unwrap();
            let spec = built.projectile.clone().unwrap();
            let max_t = built.projectile_max_time.unwrap();
            let mut sim = built.into_simulation();
            sim.inject_projectile(&spec).unwrap();
            let o = sim.run_projectile(max_t);
            let t = sim.trackers.last().unwrap();
            println!(
                "lag {lag:+5.1}: {:<12} closest gap {:>12.4e} au to {}  ({})",
                o.class().name(),
                t.closest_body.as_ref().map(|c| c.1 / AU).unwrap_or(f64::NAN),
                t.closest_body.as_ref().map(|c| c.0.as_str()).unwrap_or("-"),
                o.label()
            );
        }
    }

    #[test]
    fn figure_eight_is_periodic() {
        use crate::simulation::Simulation;
        let built = figure_eight().build().unwrap();
        let x0 = built.system.positions();
        let mut sim = Simulation::new(built.system, built.config);
        sim.advance(6.32591398);
        for (a, b) in sim.system.positions().iter().zip(&x0) {
            assert!((*a - *b).norm() < 1e-3, "{a} vs {b}");
        }
    }
}
