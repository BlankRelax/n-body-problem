//! Ready-made scenarios: the Solar System and TRAPPIST-1, Lagrange-point
//! demonstrations, periodic three-body choreographies, chaotic textbook
//! problems, rosettes, star clusters and colliding galaxies, and projectile
//! and impact experiments.

use crate::body::BodyKind;
use crate::collision::CollisionMode;
use crate::integrator::IntegratorKind;
use crate::material::Material;
use crate::orbit::{self, Elements};
use crate::projectile;
use crate::scenario::{AdaptiveSpec, BodySpec, OrbitSpec, ProjectileFileSpec, Scenario, SimulationSpec, UnitSpec, ViewSpec};
use crate::units::*;
use crate::vector::Vec3;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;

/// Every preset, with a one-line description. The order is the display order.
pub const NAMES: &[(&str, &str)] = &[
    ("solar-system", "Sun and the eight planets at J2000 (au, solar masses, years)"),
    ("inner-solar-system", "Sun, Mercury to Mars, and the Moon around Earth"),
    ("earth-moon", "Earth and Moon (km, kg, hours)"),
    ("trappist-1", "TRAPPIST-1: seven Earth-sized planets in a resonant chain, all inside 0.07 au (au, solar masses, days)"),
    ("lagrange-points", "Sun and Jupiter with probes parked at L1-L5 and Trojan swarms librating about L4 and L5; opens in the co-rotating frame"),
    ("co-orbital", "Tadpole and horseshoe orbits: test masses sharing a planet's orbit, G = 1; opens in the co-rotating frame"),
    ("lagrange-triangle", "Lagrange's equilateral solution: three unequal masses on homothetic ellipses that keep forming the same triangle, G = 1"),
    ("figure-eight", "Chenciner-Montgomery figure-eight choreography, G = 1"),
    ("butterfly", "Suvakov-Dmitrasinovic periodic three-body orbit 'butterfly I', G = 1"),
    ("moth", "Suvakov-Dmitrasinovic periodic three-body orbit 'moth I', G = 1"),
    ("dragonfly", "Suvakov-Dmitrasinovic periodic three-body orbit 'dragonfly', G = 1"),
    ("yin-yang", "Suvakov-Dmitrasinovic periodic three-body orbit 'yin-yang I', G = 1"),
    ("goggles", "Suvakov-Dmitrasinovic periodic three-body orbit 'goggles', G = 1"),
    ("pythagorean", "Burrau's Pythagorean three-body problem (3-4-5 triangle), G = 1"),
    ("sitnikov", "Sitnikov problem: a test mass bobbing through the centre of an eccentric binary along its axis, G = 1 (the view opens tilted)"),
    ("klemperer", "Klemperer rosette: five equal stars on a rotating pentagon, an exact equilibrium that eventually breaks, G = 1"),
    ("rosette", "Eight asteroids on one eccentric orbit (e = 0.7) turned by 45 degrees each: their trails draw an eight-petal rosette"),
    ("binary-planet", "Alpha-Centauri-like binary with an Earth-mass planet around the primary"),
    ("kepler-16", "Kepler-16: a 0.69 + 0.20 solar-mass eclipsing binary with the circumbinary Saturn-mass planet Kepler-16b (au, solar masses, days)"),
    ("plummer", "200-body Plummer-sphere star cluster, Barnes-Hut forces (pc, solar masses, Myr)"),
    ("galaxy-collision", "Two disc galaxies of test stars on a slow prograde encounter: bridges and tidal tails (kpc, solar masses, Myr)"),
    ("projectile-capture", "Sun + Jupiter visited by a rogue Earth-mass planet aimed at Jupiter's orbit"),
    ("projectile-flyby", "Sun + Earth visited by a fast comet"),
    ("giant-impact", "The Moon-forming impact: Mars-sized Theia grazes the proto-Earth at 45 degrees near the escape speed (km, Earth masses, hours)"),
];

/// The presets arranged by theme, for menus.
pub const GROUPS: &[(&str, &[&str])] = &[
    ("Solar System", &["solar-system", "inner-solar-system", "earth-moon", "trappist-1"]),
    ("Lagrange points", &["lagrange-points", "co-orbital", "lagrange-triangle"]),
    ("Choreographies", &["figure-eight", "butterfly", "moth", "dragonfly", "yin-yang", "goggles"]),
    ("Chaos", &["pythagorean", "sitnikov"]),
    ("Rosettes", &["klemperer", "rosette"]),
    ("Stars and galaxies", &["binary-planet", "kepler-16", "plummer", "galaxy-collision"]),
    ("Projectiles and impacts", &["projectile-capture", "projectile-flyby", "giant-impact"]),
];

pub fn names() -> &'static [(&'static str, &'static str)] {
    NAMES
}

pub fn groups() -> &'static [(&'static str, &'static [&'static str])] {
    GROUPS
}

pub fn preset(name: &str) -> Option<Scenario> {
    Some(match name.trim().to_ascii_lowercase().as_str() {
        "solar-system" | "solar_system" | "solar" => solar_system(),
        "inner-solar-system" | "inner" => inner_solar_system(),
        "earth-moon" | "earth_moon" => earth_moon(),
        "trappist-1" | "trappist" => trappist_1(),
        "lagrange-points" | "lagrange" | "trojans" => lagrange_points(),
        "co-orbital" | "coorbital" | "horseshoe" => co_orbital(),
        "lagrange-triangle" | "triangle" => lagrange_triangle(),
        "figure-eight" | "figure8" | "figure_eight" => figure_eight(),
        "butterfly" => choreography("Butterfly", 0.30689, 0.12551, 6.2356),
        "moth" => choreography("Moth", 0.46444, 0.39606, 14.8939),
        "dragonfly" => choreography("Dragonfly", 0.08058, 0.58884, 21.2710),
        "yin-yang" | "yinyang" => choreography("Yin-yang", 0.51394, 0.30474, 17.3284),
        "goggles" => choreography("Goggles", 0.08330, 0.12789, 10.4668),
        "pythagorean" | "burrau" => pythagorean(),
        "sitnikov" => sitnikov(),
        "klemperer" => klemperer(),
        "rosette" | "flower" => rosette(),
        "binary-planet" | "binary" => binary_planet(),
        "kepler-16" | "kepler16" | "tatooine" => kepler_16(),
        "plummer" | "cluster" => plummer(200, 7),
        "galaxy-collision" | "galaxies" | "mice" => galaxy_collision(),
        "projectile-capture" | "capture" => projectile_capture(),
        "projectile-flyby" | "flyby" => projectile_flyby(),
        "giant-impact" | "theia" => giant_impact(),
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

/// A body given by an explicit state.
fn state_body(name: impl Into<String>, kind: BodyKind, mass: f64, material: Material, position: Vec3, velocity: Vec3) -> BodySpec {
    BodySpec {
        name: name.into(),
        kind: Some(kind),
        mass: Some(mass),
        material: Some(material),
        position: Some(position),
        velocity: Some(velocity),
        ..Default::default()
    }
}

fn unit_body(name: &str, mass: f64, position: [f64; 3], velocity: [f64; 3]) -> BodySpec {
    state_body(name, BodyKind::Star, mass, Material::Plasma, Vec3::from(position), Vec3::from(velocity))
}

// ---- Solar System ----------------------------------------------------------------

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
        view: ViewSpec::default(),
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
        view: ViewSpec::default(),
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
        view: ViewSpec::default(),
    }
}

/// TRAPPIST-1 (Agol et al. 2021, *Planet. Sci. J.* 2, 1): masses in Earth
/// masses, radii in Earth radii, semi-major axes in au, eccentricities.
const TRAPPIST: [(&str, f64, f64, f64, f64); 7] = [
    ("TRAPPIST-1b", 1.374, 1.116, 0.01154, 0.00622),
    ("TRAPPIST-1c", 1.308, 1.097, 0.01580, 0.00654),
    ("TRAPPIST-1d", 0.388, 0.788, 0.02227, 0.00837),
    ("TRAPPIST-1e", 0.692, 0.920, 0.02925, 0.00510),
    ("TRAPPIST-1f", 1.039, 1.045, 0.03849, 0.01007),
    ("TRAPPIST-1g", 1.321, 1.129, 0.04683, 0.00208),
    ("TRAPPIST-1h", 0.326, 0.755, 0.06189, 0.00567),
];

pub fn trappist_1() -> Scenario {
    let mut bodies = vec![BodySpec {
        name: "TRAPPIST-1".into(),
        kind: Some(BodyKind::Star),
        mass: Some(0.0898),
        radius: Some(0.1192 * R_SUN / AU),
        material: Some(Material::Plasma),
        ..Default::default()
    }];
    for (k, (name, m, r, a, e)) in TRAPPIST.iter().enumerate() {
        bodies.push(BodySpec {
            name: (*name).into(),
            kind: Some(BodyKind::Planet),
            mass: Some(m * M_EARTH / M_SUN),
            radius: Some(r * R_EARTH / AU),
            material: Some(Material::Rock),
            core: Some(Material::Iron),
            orbit: Some(OrbitSpec {
                around: "TRAPPIST-1".into(),
                a: *a,
                e: *e,
                arg_periapsis: 70.0 * k as f64,
                mean_anomaly: Some(137.5 * k as f64),
                ..Default::default()
            }),
            ..Default::default()
        });
    }
    Scenario {
        name: "TRAPPIST-1".into(),
        description: "An ultracool dwarf of 0.09 solar masses with seven Earth-sized planets whose periods run from 1.5 to 18.8 days in a chain of near resonances (8:5, 5:3, 3:2, 3:2, 4:3, 3:2). The whole system would fit inside Mercury's orbit.".into(),
        units: UnitSpec::new("au", "solar", "day"),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: Some(0.002),
            duration: Some(40.0),
            ..Default::default()
        },
        bodies,
        projectile: None,
        view: ViewSpec { trail_span: Some(19.0), ..Default::default() },
    }
}

// ---- Lagrange points -------------------------------------------------------------

/// Sun, Jupiter on a circular orbit, probes at the five Lagrange points and
/// two Trojan swarms. Everything is given co-rotating velocities about the
/// Sun's initial position; the builder's barycentric shift then makes the
/// motion exactly circular about the centre of mass.
pub fn lagrange_points() -> Scenario {
    let units = UnitSpec::astronomical();
    let u = units.resolve().expect("known units");
    let g = u.g();
    let jupiter = &PLANETS[4];
    let (m_sun, m_j) = (1.0, jupiter.mass / M_SUN);
    let a_j = 5.2029;
    let n = (g * (m_sun + m_j) / (a_j * a_j * a_j)).sqrt();
    let omega = Vec3::Z * n;
    let r_j = Vec3::new(a_j, 0.0, 0.0);
    let mut jup = planet_spec(jupiter);
    jup.orbit = None;
    jup.position = Some(r_j);
    jup.velocity = Some(omega.cross(r_j));
    let mut bodies = vec![sun_spec(), jup];
    let probe = |name: String, position: Vec3, velocity: Vec3| BodySpec { radius: Some(50e3 / AU), ..state_body(name, BodyKind::Asteroid, 1e-12, Material::Rock, position, velocity) };
    // The exact points, nudged outward by a hair so the unstable ones (L1-L3)
    // show their instability within a few Jupiter years instead of waiting for round-off.
    for (k, p) in orbit::lagrange_points(Vec3::ZERO, m_sun, r_j, m_j, Vec3::Z).iter().enumerate() {
        let position = *p * (1.0 + 2e-4 / p.norm());
        bodies.push(probe(format!("L{}", k + 1), position, omega.cross(position)));
    }
    // Trojan swarms: the Greek camp about L4 (leading) and the Trojan camp about L5 (trailing).
    const GREEKS: [&str; 10] = ["Achilles", "Hektor", "Nestor", "Agamemnon", "Odysseus", "Ajax", "Diomedes", "Menelaus", "Telamon", "Antilochus"];
    const TROJANS: [&str; 10] = ["Patroclus", "Priamus", "Aeneas", "Anchises", "Troilus", "Deiphobus", "Paris", "Antenor", "Polydoros", "Glaukos"];
    let mut rng = rand::rngs::StdRng::seed_from_u64(4);
    for (camp, centre) in [(&GREEKS, 60.0f64), (&TROJANS, -60.0f64)] {
        for name in camp.iter() {
            let theta = (centre + rng.gen_range(-7.0..7.0)).to_radians();
            let r = a_j * (1.0 + rng.gen_range(-0.01..0.01));
            let position = Vec3::new(r * theta.cos(), r * theta.sin(), rng.gen_range(-0.05..0.05));
            let velocity = omega.cross(position) * (1.0 + rng.gen_range(-0.002..0.002));
            bodies.push(probe((*name).into(), position, velocity));
        }
    }
    Scenario {
        name: "Lagrange points".into(),
        description: format!(
            "Jupiter on a circular orbit about the Sun, with a probe parked at each Lagrange point and ten Trojans around L4 and L5. In the frame co-rotating with Jupiter the probes should stand still: L4 and L5 are stable (the mass ratio {:.5} is below Routh's {:.4}) and the swarms trace tadpole-shaped librations, while L1, L2 and L3 are saddle points whose probes drift off within a few Jupiter years.",
            m_j / (m_sun + m_j),
            orbit::ROUTH_MASS_RATIO
        ),
        units,
        simulation: SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: Some(0.004),
            duration: Some(300.0),
            ..Default::default()
        },
        bodies,
        projectile: None,
        view: ViewSpec { lagrange_points: Some(["Sun".into(), "Jupiter".into()]), corotating: true, trail_span: Some(60.0), tilt: None },
    }
}

/// Test masses sharing a planet's orbit, G = 1: two tadpoles about L4 and L5
/// and two horseshoes that swing from one side of the planet to the other.
pub fn co_orbital() -> Scenario {
    let (m1, m2, a): (f64, f64, f64) = (1.0, 1e-3, 1.0);
    let n = ((m1 + m2) / (a * a * a)).sqrt();
    let omega = Vec3::Z * n;
    let r2 = Vec3::new(a, 0.0, 0.0);
    let on_circle = |theta_deg: f64, r: f64| {
        let (s, c) = theta_deg.to_radians().sin_cos();
        Vec3::new(r * c, r * s, 0.0)
    };
    let test = |name: &str, position: Vec3, velocity: Vec3| state_body(name, BodyKind::Asteroid, 1e-9, Material::Ice, position, velocity);
    // Tadpoles: at rest in the rotating frame, displaced along the orbit from L4 / L5.
    let tad_a = on_circle(60.0 + 20.0, a);
    let tad_b = on_circle(-60.0 - 45.0, a);
    // Horseshoes: start opposite the planet on slightly smaller / larger circles,
    // each with its own circular speed, so they drift round towards the planet
    // in the rotating frame and are turned back when they get there.
    let hs_in = on_circle(180.0, a * (1.0 - 0.03));
    let hs_out = on_circle(180.0, a * (1.0 + 0.03));
    let circular = |p: Vec3| Vec3::Z.cross(p).normalized() * (m1 / p.norm()).sqrt();
    Scenario {
        name: "Co-orbital motion".into(),
        description: "A planet of one thousandth of its star's mass and four test masses on the same orbit, G = 1. 'Tadpole A' and 'Tadpole B' were released at rest in the rotating frame some way from L4 and L5, so they librate about those points on tadpole-shaped paths. 'Horseshoe in' and 'Horseshoe out' start opposite the planet on orbits 3% smaller and larger; they gain on the planet, are turned around by it before they can reach it, and swing back round the other way, tracing a horseshoe in the co-rotating frame.".into(),
        units: UnitSpec::default(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            dt: Some(0.005),
            g: Some(1.0),
            duration: Some(400.0),
            collisions: CollisionMode::None,
            ..Default::default()
        },
        bodies: vec![
            unit_body("Star", m1, [0.0; 3], [0.0; 3]),
            state_body("Planet", BodyKind::Planet, m2, Material::Gas, r2, omega.cross(r2)),
            test("Tadpole A", tad_a, omega.cross(tad_a)),
            test("Tadpole B", tad_b, omega.cross(tad_b)),
            test("Horseshoe in", hs_in, circular(hs_in)),
            test("Horseshoe out", hs_out, circular(hs_out)),
        ],
        projectile: None,
        view: ViewSpec { lagrange_points: Some(["Star".into(), "Planet".into()]), corotating: true, trail_span: Some(400.0), tilt: None },
    }
}

/// Lagrange's (1772) equilateral solution: each body moves on a Kepler
/// ellipse about the centre of mass with gravitational parameter
/// `G M (ρᵢ/d)³`, all with the same eccentricity and phase, so the triangle
/// stays equilateral while it turns and breathes.
pub fn lagrange_triangle() -> Scenario {
    let masses = [1.0, 0.6, 0.3];
    let e = 0.5;
    let d = 1.0; // side length at periapsis
    let total: f64 = masses.iter().sum();
    let vertices = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(d, 0.0, 0.0), Vec3::new(0.5 * d, 0.75f64.sqrt() * d, 0.0)];
    let com = vertices.iter().zip(&masses).map(|(v, m)| *v * *m).sum::<Vec3>() / total;
    let names = ["Alpha", "Beta", "Gamma"];
    let materials = [Material::Plasma, Material::Gas, Material::Ice];
    let bodies = (0..3)
        .map(|i| {
            let rel = vertices[i] - com;
            let rho = rel.norm();
            let mu_i = total * (rho / d).powi(3); // G = 1
            let speed = (mu_i * (1.0 + e) / rho).sqrt();
            let velocity = Vec3::Z.cross(rel).normalized() * speed;
            state_body(names[i], BodyKind::Star, masses[i], materials[i].clone(), rel, velocity)
        })
        .collect();
    let period = 2.0 * PI * ((d / (1.0 - e)).powi(3) / total).sqrt();
    Scenario {
        name: "Lagrange's triangle".into(),
        description: format!(
            "Masses 1, 0.6 and 0.3 at the corners of an equilateral triangle, each on an ellipse of eccentricity 0.5 about the common centre of mass, G = 1. The triangle turns and breathes but stays equilateral (period {period:.3}). For masses this similar the configuration is unstable, so round-off eventually tips it into chaos; the stable case needs one body to carry over 96% of the mass."
        ),
        units: UnitSpec::default(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            dt: Some(0.002),
            g: Some(1.0),
            duration: Some(3.0 * period),
            collisions: CollisionMode::None,
            ..Default::default()
        },
        bodies,
        projectile: None,
        view: ViewSpec { trail_span: Some(period), ..Default::default() },
    }
}

// ---- choreographies and chaos -------------------------------------------------------

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
        view: ViewSpec { trail_span: Some(6.32591398), ..Default::default() },
    }
}

/// One of the planar periodic three-body orbits of Šuvakov & Dmitrašinović
/// (2013, *PRL* 110, 114301) in their convention: equal masses, G = 1, bodies
/// at (−1, 0), (1, 0) and (0, 0), the outer two with velocity `(vx, vy)` and
/// the middle one with `(−2vx, −2vy)`. The published initial conditions have
/// five digits, so the orbits are periodic to about that precision and drift
/// after several periods.
pub fn choreography(title: &str, vx: f64, vy: f64, period: f64) -> Scenario {
    Scenario {
        name: format!("Three-body '{title}'"),
        description: format!(
            "Šuvakov & Dmitrašinović's periodic orbit '{}' (period {period:.4}): three equal masses, G = 1. Watch one full period draw itself, then another; the published five-digit initial conditions let it wander after several periods.",
            title.to_lowercase()
        ),
        units: UnitSpec::default(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            adaptive: Some(AdaptiveSpec { eta: 0.005, dt_min: 1e-9, dt_max: 0.004 }),
            g: Some(1.0),
            duration: Some(2.0 * period),
            collisions: CollisionMode::None,
            ..Default::default()
        },
        bodies: vec![
            unit_body("A", 1.0, [-1.0, 0.0, 0.0], [vx, vy, 0.0]),
            unit_body("B", 1.0, [1.0, 0.0, 0.0], [vx, vy, 0.0]),
            unit_body("C", 1.0, [0.0, 0.0, 0.0], [-2.0 * vx, -2.0 * vy, 0.0]),
        ],
        projectile: None,
        view: ViewSpec { trail_span: Some(period), ..Default::default() },
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
        view: ViewSpec::default(),
    }
}

/// Sitnikov (1960): a massless body on the axis of symmetry of an eccentric
/// equal-mass binary. Its bobbing is chaotic for large amplitudes.
pub fn sitnikov() -> Scenario {
    let (m, a, e): (f64, f64, f64) = (0.5, 1.0, 0.6);
    // Relative orbit at apoapsis: separation a(1+e), speed sqrt(μ(1−e)/(a(1+e))) with μ = 2m.
    let sep = a * (1.0 + e);
    let v_rel = (2.0 * m * (1.0 - e) / sep).sqrt();
    Scenario {
        name: "Sitnikov problem".into(),
        description: "Two equal stars on an orbit of eccentricity 0.6, and a test mass dropped from rest on the axis through their centre of mass, G = 1. The test mass falls through the plane and back, but the stars' changing separation pumps or drains its energy each time, so its bobbing is chaotic and it can be thrown off altogether. Tilt the view (right-drag) to see the third dimension.".into(),
        units: UnitSpec::default(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            adaptive: Some(AdaptiveSpec { eta: 0.01, dt_min: 1e-9, dt_max: 0.005 }),
            g: Some(1.0),
            duration: Some(80.0),
            collisions: CollisionMode::None,
            ..Default::default()
        },
        bodies: vec![
            unit_body("Star A", m, [-0.5 * sep, 0.0, 0.0], [0.0, -0.5 * v_rel, 0.0]),
            unit_body("Star B", m, [0.5 * sep, 0.0, 0.0], [0.0, 0.5 * v_rel, 0.0]),
            state_body("Bob", BodyKind::Asteroid, 1e-9, Material::Ice, Vec3::new(0.0, 0.0, 1.3), Vec3::ZERO),
        ],
        projectile: None,
        view: ViewSpec { tilt: Some(65.0), trail_span: Some(25.0), ..Default::default() },
    }
}

// ---- rosettes -----------------------------------------------------------------------

/// Klemperer (1962): `n` equal masses at the corners of a regular polygon,
/// rotating rigidly. The inward pull on each body is
/// `(G m / 4R²) Σₖ 1/sin(πk/n)`, which fixes the circular speed.
pub fn klemperer() -> Scenario {
    let (n, m, r) = (5usize, 1.0, 1.0);
    let pull: f64 = (1..n).map(|k| 1.0 / (PI * k as f64 / n as f64).sin()).sum::<f64>() * m / (4.0 * r * r);
    let speed = (pull * r).sqrt();
    let period = 2.0 * PI * r / speed;
    let materials = [Material::Plasma, Material::Ice, Material::Gas, Material::Water, Material::Rubber];
    let bodies = (0..n)
        .map(|k| {
            let theta = 2.0 * PI * k as f64 / n as f64;
            let position = Vec3::new(r * theta.cos(), r * theta.sin(), 0.0);
            state_body(format!("Star {}", k + 1), BodyKind::Star, m, materials[k % materials.len()].clone(), position, Vec3::Z.cross(position).normalized() * speed)
        })
        .collect();
    Scenario {
        name: "Klemperer rosette".into(),
        description: format!(
            "Five equal stars at the corners of a regular pentagon, each moving at exactly the speed that makes the pull of the other four hold it on its circle (period {period:.3}), G = 1. The configuration is an exact solution but an unstable one: rounding errors of one part in 10^16 grow exponentially until, after several turns, the pentagon tears itself apart."
        ),
        units: UnitSpec::default(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            adaptive: Some(AdaptiveSpec { eta: 0.005, dt_min: 1e-9, dt_max: 0.005 }),
            g: Some(1.0),
            duration: Some(20.0 * period),
            collisions: CollisionMode::Merge,
            ..Default::default()
        },
        bodies,
        projectile: None,
        view: ViewSpec { trail_span: Some(2.0 * period), ..Default::default() },
    }
}

/// Eight asteroids on the same eccentric orbit about a Sun-like star, each
/// orbit turned by 45° from the last.
pub fn rosette() -> Scenario {
    let n = 8;
    let mut bodies = vec![sun_spec()];
    for k in 0..n {
        let turn = 360.0 * k as f64 / n as f64;
        bodies.push(BodySpec {
            name: format!("Petal {}", k + 1),
            kind: Some(BodyKind::Asteroid),
            mass: Some(1e-10),
            material: Some(if k % 2 == 0 { Material::Ice } else { Material::Rock }),
            orbit: Some(OrbitSpec { around: "Sun".into(), a: 1.0, e: 0.7, arg_periapsis: turn, mean_anomaly: Some(turn * 3.0), ..Default::default() }),
            ..Default::default()
        });
    }
    Scenario {
        name: "Rosette".into(),
        description: "Eight Ceres-mass asteroids on one-year orbits of eccentricity 0.7 about the Sun, each orbit turned by 45 degrees from its neighbour. Over a year their trails draw an eight-petal rosette; the asteroids speed through periapsis at 0.3 au and dawdle at 1.7 au.".into(),
        units: UnitSpec::astronomical(),
        simulation: SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: Some(0.0005),
            duration: Some(2.0),
            ..Default::default()
        },
        bodies,
        projectile: None,
        view: ViewSpec { trail_span: Some(1.0), ..Default::default() },
    }
}

// ---- stars and galaxies -------------------------------------------------------------

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
        view: ViewSpec::default(),
    }
}

/// Kepler-16 (Doyle et al. 2011, *Science* 333, 1602): a circumbinary planet.
/// The planet orbits the pair's centre of mass, so its state is computed
/// explicitly rather than relative to one star.
pub fn kepler_16() -> Scenario {
    let units = UnitSpec::new("au", "solar", "day");
    let u = units.resolve().expect("known units");
    let g = u.g();
    let (m_a, m_b, m_p) = (0.6897, 0.20255, 0.333 * M_JUPITER / M_SUN);
    // Binary: a = 0.22431 au, e = 0.15944; planet: a = 0.7048 au, e = 0.0069.
    let el_b = Elements { a: 0.22431, e: 0.15944, inclination: 0.0, raan: 0.0, arg_periapsis: 0.0, true_anomaly: 0.5 * PI };
    let (r_b, v_b) = orbit::state_from_elements(g * (m_a + m_b), &el_b).expect("valid elements");
    let bary = r_b * (m_b / (m_a + m_b));
    let v_bary = v_b * (m_b / (m_a + m_b));
    let el_p = Elements { a: 0.7048, e: 0.0069, inclination: 0.0, raan: 0.0, arg_periapsis: 0.0, true_anomaly: PI };
    let (r_p, v_p) = orbit::state_from_elements(g * (m_a + m_b + m_p), &el_p).expect("valid elements");
    Scenario {
        name: "Kepler-16".into(),
        description: "The first circumbinary planet seen to transit: a K dwarf and an M dwarf 0.22 au apart (period 41 days) circled by a Saturn-mass planet at 0.70 au (229 days). Both stars orbit their common centre of mass, so the planet sees its suns swing past each other every three weeks.".into(),
        units,
        simulation: SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: Some(0.05),
            duration: Some(700.0),
            ..Default::default()
        },
        bodies: vec![
            BodySpec { radius: Some(0.6489 * R_SUN / AU), ..state_body("Kepler-16 A", BodyKind::Star, m_a, Material::Plasma, Vec3::ZERO, Vec3::ZERO) },
            BodySpec { radius: Some(0.22623 * R_SUN / AU), ..state_body("Kepler-16 B", BodyKind::Star, m_b, Material::Plasma, r_b, v_b) },
            BodySpec {
                radius: Some(0.7538 * R_JUPITER / AU),
                core: Some(Material::Rock),
                ..state_body("Kepler-16b", BodyKind::Planet, m_p, Material::Gas, bary + r_p, v_bary + v_p)
            },
        ],
        projectile: None,
        view: ViewSpec { trail_span: Some(229.0), ..Default::default() },
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
        let phi = 2.0 * PI * rng.gen::<f64>();
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
        bodies.push(state_body(format!("star{:03}", k + 1), BodyKind::Star, total_mass / n as f64, Material::Plasma, position, velocity));
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
        view: ViewSpec::default(),
    }
}

/// Toomre & Toomre (1972) style: two point-mass nuclei, each dressed with
/// rings of test stars on circular orbits, on a slow prograde encounter.
pub fn galaxy_collision() -> Scenario {
    let units = UnitSpec::new("kpc", "solar", "myr");
    let u = units.resolve().expect("known units");
    let g = u.g();
    let (m1, m2) = (1.0e10, 0.5e10);
    let star_mass = 1.0e4;
    // Relative orbit of the nuclei: e = 0.9 (slightly bound, so they merge in the end), periapsis 6 kpc, starting 25 kpc apart on the way in.
    let (e, q, r0): (f64, f64, f64) = (0.9, 6.0, 25.0);
    let a = q / (1.0 - e);
    let p = a * (1.0 - e * e);
    let nu0 = -((p / r0 - 1.0) / e).clamp(-1.0, 1.0).acos();
    let el = Elements { a, e, inclination: 0.0, raan: 0.0, arg_periapsis: 0.0, true_anomaly: nu0 };
    let (r_rel, v_rel) = orbit::state_from_elements(g * (m1 + m2), &el).expect("valid elements");
    let mut bodies = vec![
        state_body("Nucleus A", BodyKind::Star, m1, Material::Plasma, Vec3::ZERO, Vec3::ZERO),
        state_body("Nucleus B", BodyKind::Star, m2, Material::Plasma, r_rel, v_rel),
    ];
    let mut disc = |label: &str, centre: Vec3, v_centre: Vec3, mass: f64, rings: usize, r_in: f64, r_out: f64, tilt_deg: f64| {
        let (st, ct) = tilt_deg.to_radians().sin_cos();
        let tilt = |v: Vec3| Vec3::new(v.x, ct * v.y - st * v.z, st * v.y + ct * v.z);
        let mut k = 0;
        for ring in 0..rings {
            let r = r_in + (r_out - r_in) * ring as f64 / (rings - 1).max(1) as f64;
            let count = 8 + 3 * ring;
            let speed = (g * mass / r).sqrt();
            for j in 0..count {
                let phi = 2.0 * PI * j as f64 / count as f64 + 0.37 * ring as f64;
                let pos = Vec3::new(r * phi.cos(), r * phi.sin(), 0.0);
                let vel = Vec3::Z.cross(pos).normalized() * speed;
                k += 1;
                bodies.push(state_body(format!("{label}{k:03}"), BodyKind::Particle, star_mass, if label == "a" { Material::Ice } else { Material::Gas }, centre + tilt(pos), v_centre + tilt(vel)));
            }
        }
    };
    disc("a", Vec3::ZERO, Vec3::ZERO, m1, 8, 1.5, 7.0, 0.0);
    disc("b", r_rel, v_rel, m2, 6, 1.2, 5.0, 30.0);
    Scenario {
        name: "Galaxy collision".into(),
        description: "Two disc galaxies, 10 and 5 billion solar masses, each a point-mass nucleus dressed with rings of test stars (148 and 93) on circular orbits, approaching on a slow prograde orbit with a closest approach of 6 kpc about 300 Myr in. The passage draws out tidal tails and a bridge between them, in the manner of Toomre & Toomre (1972); being slightly bound, the pair falls back together and merges.".into(),
        units,
        simulation: SimulationSpec {
            integrator: IntegratorKind::Yoshida4,
            dt: Some(0.5),
            softening: 0.05,
            force: "barnes_hut".into(),
            theta: 0.6,
            collisions: CollisionMode::None,
            duration: Some(1500.0),
            ..Default::default()
        },
        bodies,
        projectile: None,
        view: ViewSpec { trail_span: Some(250.0), ..Default::default() },
    }
}

// ---- projectiles and impacts -------------------------------------------------------

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
    let nu_j = (phi - n_j * geo.time_to_periapsis + lag_deg.to_radians()).rem_euclid(2.0 * PI).to_degrees();
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
        view: ViewSpec::default(),
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
        view: ViewSpec::default(),
    }
}

/// The canonical Moon-forming impact (Canup 2004): a Mars-mass impactor
/// hitting the proto-Earth obliquely at close to the mutual escape speed.
/// Here it is set up as a projectile so the collision model can be watched
/// deciding between merging and a hit-and-run.
pub fn giant_impact() -> Scenario {
    let units = UnitSpec::new("km", "earth", "hour");
    let u = units.resolve().expect("known units");
    let g = u.g();
    let (m_earth, r_earth) = (0.9, 6250.0);
    let (m_theia, r_theia) = (0.1, 3300.0);
    let mu = g * (m_earth + m_theia);
    let r_contact = r_earth + r_theia;
    let v_esc = (2.0 * mu / r_contact).sqrt();
    let v_inf = 0.3 * v_esc;
    let v_contact = (v_inf * v_inf + v_esc * v_esc).sqrt();
    let angle = 45f64.to_radians();
    // Angular momentum: b v∞ = r v sin θ at contact.
    let b = r_contact * v_contact * angle.sin() / v_inf;
    Scenario {
        name: "Giant impact".into(),
        description: format!(
            "A Mars-sized body of a tenth of an Earth mass strikes the 0.9-Earth-mass proto-Earth at a 45-degree angle and {:.2} times their mutual escape speed ({:.1} km/s at contact), the canonical Moon-forming impact. Press Launch: the physical collision model decides whether the pair merges or Theia survives a hit-and-run; the merger's layers differentiate, iron sinking to the core.",
            v_contact / v_esc,
            v_contact * u.velocity() / 1e3
        ),
        units,
        simulation: SimulationSpec {
            integrator: IntegratorKind::Pefrl,
            adaptive: Some(AdaptiveSpec { eta: 0.01, dt_min: 1e-7, dt_max: 0.01 }),
            duration: Some(40.0),
            ..Default::default()
        },
        bodies: vec![BodySpec {
            name: "Proto-Earth".into(),
            kind: Some(BodyKind::Planet),
            mass: Some(m_earth),
            radius: Some(r_earth),
            material: Some(Material::Rock),
            core: Some(Material::Iron),
            ..Default::default()
        }],
        projectile: Some(ProjectileFileSpec {
            name: "Theia".into(),
            mass: Some(m_theia),
            radius: Some(r_theia),
            material: Some(Material::Rock),
            core: Some(Material::Iron),
            v_infinity: v_inf,
            impact_parameter: b,
            direction: Vec3::X,
            start_distance: Some(12.0 * r_contact),
            max_time: Some(40.0),
            ..Default::default()
        }),
        view: ViewSpec { trail_span: Some(12.0), ..Default::default() },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::Diagnostics;
    use crate::simulation::Simulation;

    #[test]
    fn all_presets_build() {
        for (name, _) in NAMES {
            let sc = preset(name).unwrap_or_else(|| panic!("missing preset {name}"));
            let built = sc.build().unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(!built.system.is_empty());
            let text = sc.to_toml().unwrap();
            let again = Scenario::from_toml(&text).unwrap();
            assert_eq!(again.view, sc.view, "{name}: view hints round trip");
            again.build().unwrap_or_else(|e| panic!("{name} round trip: {e}"));
        }
    }

    #[test]
    fn groups_cover_every_preset_once() {
        let mut seen = Vec::new();
        for (_, members) in GROUPS {
            for m in members.iter() {
                assert!(NAMES.iter().any(|(n, _)| n == m), "group member {m} is not a preset");
                assert!(!seen.contains(m), "{m} listed twice");
                seen.push(*m);
            }
        }
        for (n, _) in NAMES {
            assert!(seen.contains(n), "{n} is in no group");
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
    fn trappist_periods_match_the_discovery_paper() {
        let built = trappist_1().build().unwrap();
        let star = built.system.find("TRAPPIST-1").unwrap();
        let periods = [1.51088, 2.42179, 4.04978, 6.09956, 9.20659, 12.3535, 18.7672];
        for ((name, ..), expected) in TRAPPIST.iter().zip(periods) {
            let b = built.system.find(name).unwrap();
            let mu = built.system.g * (star.mass() + b.mass());
            let el = crate::orbit::elements_from_state(mu, b.position - star.position, b.velocity - star.velocity);
            let days = el.period.unwrap() / DAY;
            assert!((days / expected - 1.0).abs() < 4e-3, "{name}: {days} d vs {expected} d");
        }
    }

    #[test]
    fn plummer_is_near_virial_equilibrium() {
        let built = plummer(400, 3).build().unwrap();
        let d = Diagnostics::compute(&built.system);
        assert!((d.virial_ratio - 1.0).abs() < 0.25, "virial ratio {}", d.virial_ratio);
        assert!(d.total < 0.0);
    }

    /// Longitude of `body` measured from `secondary` in the frame co-rotating
    /// with the `primary`-`secondary` pair, degrees in (−180, 180].
    fn corotating_longitude(system: &crate::system::NBodySystem, primary: &str, secondary: &str, body: &str) -> f64 {
        let p = system.find(primary).unwrap().position;
        let s = system.find(secondary).unwrap().position - p;
        let b = system.find(body).unwrap().position - p;
        let d = b.y.atan2(b.x) - s.y.atan2(s.x);
        (d.to_degrees() + 180.0).rem_euclid(360.0) - 180.0
    }

    /// Probes at L4 and L5 stay put, the Trojans stay in their camps, and the
    /// collinear probes leave within a few Jupiter years.
    #[test]
    fn lagrange_probes_behave() {
        let built = lagrange_points().build().unwrap();
        let sun_pos = built.system.find("Sun").unwrap().position;
        let a_j = built.system.find("Jupiter").unwrap().position.distance(sun_pos);
        let hill = a_j * (built_jupiter_mass() / 3.0f64).cbrt();
        let mut sim = Simulation::new(built.system, built.config);
        let mut max_dev = [0.0f64; 2];
        let mut camp_dev = 0.0f64;
        let mut collinear_drift = [0.0f64; 3];
        for _ in 0..60 {
            sim.advance(5.0 * YEAR);
            let (sun, jup) = (sim.system.find("Sun").unwrap().clone(), sim.system.find("Jupiter").unwrap().clone());
            let pts = crate::orbit::lagrange_points(sun.position, sun.mass(), jup.position, jup.mass(), Vec3::Z);
            for k in 0..3 {
                let probe = sim.system.find(&format!("L{}", k + 1)).unwrap().position;
                collinear_drift[k] = collinear_drift[k].max(probe.distance(pts[k]));
            }
            for (k, name) in ["L4", "L5"].iter().enumerate() {
                let lon = corotating_longitude(&sim.system, "Sun", "Jupiter", name);
                let expected = if k == 0 { 60.0 } else { -60.0 };
                max_dev[k] = max_dev[k].max((lon - expected).abs());
            }
            for b in sim.system.bodies() {
                if b.name.len() <= 2 || b.name == "Sun" || b.name == "Jupiter" {
                    continue;
                }
                let lon = corotating_longitude(&sim.system, "Sun", "Jupiter", &b.name);
                let centre = if lon > 0.0 { 60.0 } else { -60.0 };
                camp_dev = camp_dev.max((lon - centre).abs());
                assert!(b.position.distance(jup.position) > 0.3 * a_j, "{} came within {} au of Jupiter", b.name, b.position.distance(jup.position) / AU);
            }
        }
        assert!(max_dev[0] < 2.0 && max_dev[1] < 2.0, "L4/L5 probes wandered by {max_dev:?} degrees over 300 yr");
        // Tadpoles of up to ±70° about L4/L5 (real Trojans librate by up to ~35°); none crosses over to the other camp.
        assert!(camp_dev < 100.0, "a Trojan wandered {camp_dev} degrees from its camp");
        // The collinear points are saddles: their probes leave within the run (L3 is the slowest, with an e-folding time of decades).
        for (k, d) in collinear_drift.iter().enumerate() {
            assert!(*d > 0.5 * hill, "L{} probe only drifted {} Hill radii in 300 yr", k + 1, d / hill);
        }
    }

    fn built_jupiter_mass() -> f64 {
        PLANETS[4].mass / M_SUN
    }

    /// Tadpoles librate about L4/L5; horseshoes reverse their drift before
    /// reaching the planet and visit both sides of it.
    #[test]
    fn co_orbital_tadpoles_and_horseshoes() {
        let built = co_orbital().build().unwrap();
        let mut sim = Simulation::new(built.system, built.config);
        let mut tad = [(f64::INFINITY, f64::NEG_INFINITY); 2];
        let mut hs = [(f64::INFINITY, f64::NEG_INFINITY); 2];
        let mut closest = f64::INFINITY;
        for _ in 0..400 {
            sim.advance(1.0);
            for (k, name) in ["Tadpole A", "Tadpole B"].iter().enumerate() {
                let lon = corotating_longitude(&sim.system, "Star", "Planet", name);
                tad[k] = (tad[k].0.min(lon), tad[k].1.max(lon));
            }
            for (k, name) in ["Horseshoe in", "Horseshoe out"].iter().enumerate() {
                let lon = corotating_longitude(&sim.system, "Star", "Planet", name);
                hs[k] = (hs[k].0.min(lon), hs[k].1.max(lon));
                closest = closest.min(sim.system.find(name).unwrap().position.distance(sim.system.find("Planet").unwrap().position));
            }
        }
        assert!(tad[0].0 > 20.0 && tad[0].1 < 180.0, "tadpole A longitude range {:?}", tad[0]);
        assert!(tad[1].1 < -20.0 && tad[1].0 > -180.0, "tadpole B longitude range {:?}", tad[1]);
        for (k, range) in hs.iter().enumerate() {
            assert!(range.0 < -90.0 && range.1 > 90.0, "horseshoe {k} did not visit both sides: {range:?}");
        }
        assert!(closest > 0.15, "a horseshoe came within {closest} of the planet");
    }

    #[test]
    fn lagrange_triangle_stays_equilateral() {
        let built = lagrange_triangle().build().unwrap();
        let mut sim = Simulation::new(built.system, built.config);
        let period = 2.0 * PI * ((1.0f64 / 0.5).powi(3) / 1.9).sqrt();
        let x0 = sim.system.positions();
        for _ in 0..40 {
            sim.advance(period / 20.0);
            let p = sim.system.positions();
            let sides = [p[0].distance(p[1]), p[1].distance(p[2]), p[2].distance(p[0])];
            let mean = sides.iter().sum::<f64>() / 3.0;
            for s in sides {
                assert!((s / mean - 1.0).abs() < 1e-5, "sides {sides:?} at t = {}", sim.system.time);
            }
        }
        for (a, b) in sim.system.positions().iter().zip(&x0) {
            assert!((*a - *b).norm() < 1e-3, "not periodic: {a} vs {b}");
        }
    }

    #[test]
    fn choreographies_close_after_one_period() {
        for (name, period) in [("butterfly", 6.2356), ("moth", 14.8939), ("dragonfly", 21.2710), ("yin-yang", 17.3284), ("goggles", 10.4668)] {
            let built = preset(name).unwrap().build().unwrap();
            let x0 = built.system.positions();
            let mut sim = Simulation::new(built.system, built.config);
            sim.advance(period);
            let err = sim.system.positions().iter().zip(&x0).map(|(a, b)| (*a - *b).norm()).fold(0.0, f64::max);
            assert!(err < 0.05, "{name}: returns to within {err} after one period");
        }
    }

    #[test]
    fn klemperer_rotates_rigidly_before_breaking() {
        let built = klemperer().build().unwrap();
        let mut sim = Simulation::new(built.system, built.config);
        let r0 = sim.system.positions()[0].norm();
        for _ in 0..20 {
            sim.advance(0.5);
            let p = sim.system.positions();
            for k in 0..5 {
                assert!((p[k].norm() / r0 - 1.0).abs() < 1e-3, "radius drifted at t = {}", sim.system.time);
                let side = p[k].distance(p[(k + 1) % 5]);
                assert!((side / (2.0 * r0 * (PI / 5.0).sin()) - 1.0).abs() < 1e-3, "side drifted");
            }
        }
    }

    /// Print each Trojan's libration range. `cargo test -p nbody-core trojan_ranges -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn trojan_ranges() {
        let built = lagrange_points().build().unwrap();
        let names: Vec<String> = built.system.bodies().iter().map(|b| b.name.clone()).filter(|n| n.len() > 2 && n != "Sun" && n != "Jupiter").collect();
        let mut sim = Simulation::new(built.system, built.config);
        let mut range: Vec<(f64, f64)> = vec![(f64::INFINITY, f64::NEG_INFINITY); names.len()];
        for _ in 0..120 {
            sim.advance(2.5 * YEAR);
            for (k, n) in names.iter().enumerate() {
                let lon = corotating_longitude(&sim.system, "Sun", "Jupiter", n);
                range[k] = (range[k].0.min(lon), range[k].1.max(lon));
            }
        }
        for (n, r) in names.iter().zip(&range) {
            println!("{n:<12} {:+7.1} .. {:+7.1}", r.0, r.1);
        }
    }

    /// Print how long the rosette survives. `cargo test -p nbody-core klemperer_lifetime -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn klemperer_lifetime() {
        let built = klemperer().build().unwrap();
        let mut sim = Simulation::new(built.system, built.config);
        let r0 = sim.system.positions()[0].norm();
        loop {
            sim.advance(1.0);
            let p = sim.system.positions();
            let dev = (0..5).map(|k| (p[k].norm() / r0 - 1.0).abs()).fold(0.0, f64::max);
            println!("t = {:6.1}  radius deviation {dev:.3e}", sim.system.time);
            if dev > 0.1 || sim.system.time > 200.0 {
                break;
            }
        }
    }

    #[test]
    fn kepler_16_planet_orbits_the_barycentre() {
        let built = kepler_16().build().unwrap();
        let mut sim = Simulation::new(built.system, built.config);
        let a_expected = 0.7048 * AU;
        for _ in 0..10 {
            sim.advance(23.0 * DAY);
            let s = &sim.system;
            let (a, b, p) = (s.find("Kepler-16 A").unwrap(), s.find("Kepler-16 B").unwrap(), s.find("Kepler-16b").unwrap());
            let m = a.mass() + b.mass();
            let bary = (a.position * a.mass() + b.position * b.mass()) / m;
            let vbary = (a.velocity * a.mass() + b.velocity * b.mass()) / m;
            let el = crate::orbit::elements_from_state(s.g * (m + p.mass()), p.position - bary, p.velocity - vbary);
            assert!((el.elements.a / a_expected - 1.0).abs() < 0.02, "a = {} au", el.elements.a / AU);
            // The binary forces an osculating eccentricity of a few per cent on top of the planet's own 0.0069.
            assert!(el.elements.e < 0.1, "e = {}", el.elements.e);
        }
    }

    #[test]
    fn galaxy_discs_start_bound_to_their_nuclei() {
        let built = galaxy_collision().build().unwrap();
        let s = &built.system;
        for (label, nucleus) in [("a", "Nucleus A"), ("b", "Nucleus B")] {
            let n = s.find(nucleus).unwrap();
            for b in s.bodies().iter().filter(|b| b.name.starts_with(label) && b.name.len() == 4) {
                let r = b.position - n.position;
                let v = b.velocity - n.velocity;
                let eps = 0.5 * v.norm_sq() - s.g * n.mass() / r.norm();
                assert!(eps < 0.0, "{} is not bound to {nucleus}", b.name);
            }
        }
        assert_eq!(s.len(), 2 + 148 + 93);
    }

    #[test]
    fn giant_impact_collides() {
        let built = giant_impact().build().unwrap();
        let spec = built.projectile.clone().unwrap();
        let max_t = built.projectile_max_time.unwrap();
        let mut sim = built.into_simulation();
        sim.inject_projectile(&spec).unwrap();
        sim.run_projectile(max_t);
        assert!(!sim.collisions.is_empty(), "Theia missed");
        let e = &sim.collisions[0];
        assert!((e.impact_angle.to_degrees() - 45.0).abs() < 5.0, "impact angle {}", e.impact_angle.to_degrees());
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
        let built = figure_eight().build().unwrap();
        let x0 = built.system.positions();
        let mut sim = Simulation::new(built.system, built.config);
        sim.advance(6.32591398);
        for (a, b) in sim.system.positions().iter().zip(&x0) {
            assert!((*a - *b).norm() < 1e-3, "{a} vs {b}");
        }
    }
}
