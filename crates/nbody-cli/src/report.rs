//! Human-readable tables and unit-aware formatting.

use nbody_core::body::Body;
use nbody_core::orbit;
use nbody_core::projectile::{Outcome, Phase, ProjectileTracker};
use nbody_core::scenario::UnitSpec;
use nbody_core::simulation::{Simulation, TimestepPolicy};
use nbody_core::system::NBodySystem;
use nbody_core::units::UnitSystem;
use nbody_core::vector::Vec3;
use std::fmt::Write as _;

/// Compact number: 6 significant digits, scientific outside [1e-3, 1e6).
pub fn sig(x: f64) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    if !x.is_finite() {
        return format!("{x}");
    }
    let a = x.abs();
    if (1e-3..1e6).contains(&a) {
        let decimals = (5 - a.log10().floor() as i32).max(0) as usize;
        let s = format!("{x:.decimals$}");
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    } else {
        format!("{x:.4e}")
    }
}

#[derive(Debug, Clone)]
pub struct Units {
    pub sys: UnitSystem,
    pub spec: UnitSpec,
}

impl Units {
    pub fn new(spec: UnitSpec) -> anyhow::Result<Units> {
        let sys = spec.resolve()?;
        Ok(Units { sys, spec })
    }

    pub fn len(&self, m: f64) -> String {
        format!("{} {}", sig(m / self.sys.length), self.spec.length)
    }

    pub fn mass(&self, kg: f64) -> String {
        format!("{} {}", sig(kg / self.sys.mass), self.spec.mass)
    }

    pub fn time(&self, s: f64) -> String {
        format!("{} {}", sig(s / self.sys.time), self.spec.time)
    }

    pub fn speed(&self, v: f64) -> String {
        format!("{} {}/{}", sig(v / self.sys.velocity()), self.spec.length, self.spec.time)
    }

    pub fn vec_len(&self, v: Vec3) -> String {
        let s = v / self.sys.length;
        format!("[{}, {}, {}] {}", sig(s.x), sig(s.y), sig(s.z), self.spec.length)
    }

    pub fn vec_speed(&self, v: Vec3) -> String {
        let s = v / self.sys.velocity();
        format!("[{}, {}, {}] {}/{}", sig(s.x), sig(s.y), sig(s.z), self.spec.length, self.spec.time)
    }
}

pub fn body_table(system: &NBodySystem, u: &Units) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "{:<3} {:<14} {:<10} {:>14} {:>14} {:<12} {:>5}  {:<38} {}",
        "#",
        "name",
        "kind",
        format!("mass [{}]", u.spec.mass),
        format!("radius [{}]", u.spec.length),
        "material",
        "shell",
        format!("position [{}]", u.spec.length),
        format!("velocity [{}/{}]", u.spec.length, u.spec.time)
    );
    for (i, b) in system.bodies().iter().enumerate() {
        let p = b.position / u.sys.length;
        let v = b.velocity / u.sys.velocity();
        let _ = writeln!(
            out,
            "{:<3} {:<14} {:<10} {:>14} {:>14} {:<12} {:>5}  {:<38} {}",
            i + 1,
            truncate(&b.name, 14),
            b.kind.name(),
            sig(b.mass() / u.sys.mass),
            sig(b.radius() / u.sys.length),
            truncate(b.material().name(), 12),
            b.structure().layers().len(),
            format!("[{}, {}, {}]", sig(p.x), sig(p.y), sig(p.z)),
            format!("[{}, {}, {}]", sig(v.x), sig(v.y), sig(v.z)),
        );
    }
    let _ = writeln!(
        out,
        "{} bodies, total mass {}, t = {}, G = {:.6e} SI",
        system.len(),
        u.mass(system.total_mass()),
        u.time(system.time),
        system.g
    );
    out
}

fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() <= n {
        s.to_string()
    } else {
        s.chars().take(n - 1).collect::<String>() + "…"
    }
}

pub fn body_detail(system: &NBodySystem, body: &Body, u: &Units) -> String {
    let mut out = String::new();
    let s = body.structure();
    let _ = writeln!(out, "{} ({}, id {})", body.name, body.kind.name(), body.id());
    let _ = writeln!(out, "  mass            {}  ({:.4e} kg)", u.mass(body.mass()), body.mass());
    let _ = writeln!(out, "  radius          {}  ({:.4e} m)", u.len(body.radius()), body.radius());
    let _ = writeln!(out, "  mean density    {:.1} kg/m^3", s.mean_density());
    let _ = writeln!(out, "  material        {} (dominant), restitution {:.2}", body.material().name(), body.restitution());
    let _ = writeln!(out, "  I / (M R^2)     {:.4}  (0.4 = uniform sphere)", s.inertia_factor());
    let _ = writeln!(out, "  escape speed    {}", u.speed(body.escape_velocity(system.g)));
    let _ = writeln!(out, "  surface gravity {:.4} m/s^2", body.surface_gravity(system.g));
    let _ = writeln!(out, "  layers (inside out):");
    let total = s.mass();
    let mut inner = 0.0;
    for l in s.layers() {
        let m = l.material.density() * (nbody_core::structure::sphere_volume(l.outer_radius) - nbody_core::structure::sphere_volume(inner));
        let _ = writeln!(
            out,
            "    r = {:>14}  {:<10} rho = {:>7.0} kg/m^3  {:>5.1}% of mass",
            u.len(l.outer_radius),
            l.material.name(),
            l.material.density(),
            100.0 * m / total
        );
        inner = l.outer_radius;
    }
    let _ = writeln!(out, "  position        {}", u.vec_len(body.position));
    let _ = writeln!(out, "  velocity        {}  (|v| = {})", u.vec_speed(body.velocity), u.speed(body.velocity.norm()));
    if let Some(primary) = system.most_massive() {
        if primary.id() != body.id() {
            let mu = system.g * (primary.mass() + body.mass());
            let r = body.position - primary.position;
            let v = body.velocity - primary.velocity;
            let el = orbit::elements_from_state(mu, r, v);
            let _ = writeln!(out, "  relative to {} (the most massive body):", primary.name);
            let _ = writeln!(out, "    distance      {}", u.len(r.norm()));
            let _ = writeln!(out, "    {}", elements_line(&el, u));
            let hill = orbit::hill_radius(el.elements.a.abs(), el.elements.e.min(0.999), body.mass(), primary.mass());
            let _ = writeln!(out, "    Hill radius   {}", u.len(hill));
            let roche_r = orbit::roche_limit_rigid(primary.radius(), primary.mean_density(), body.mean_density());
            let roche_f = orbit::roche_limit_fluid(primary.radius(), primary.mean_density(), body.mean_density());
            let flag = if r.norm() < roche_r { "  <-- INSIDE rigid Roche limit" } else if r.norm() < roche_f { "  <-- inside fluid Roche limit" } else { "" };
            let _ = writeln!(out, "    Roche limit   rigid {}, fluid {}{}", u.len(roche_r), u.len(roche_f), flag);
        }
    }
    out
}

pub fn elements_line(el: &orbit::OrbitalElements, u: &Units) -> String {
    let e = &el.elements;
    match el.kind {
        orbit::OrbitKind::Elliptic => format!(
            "a = {}, e = {:.5}, i = {:.3}°, Ω = {:.2}°, ω = {:.2}°, ν = {:.2}°, T = {}, q = {}, Q = {}",
            u.len(e.a),
            e.e,
            e.inclination.to_degrees(),
            e.raan.to_degrees(),
            e.arg_periapsis.to_degrees(),
            e.true_anomaly.to_degrees(),
            u.time(el.period.unwrap_or(f64::NAN)),
            u.len(el.periapsis),
            u.len(el.apoapsis.unwrap_or(f64::NAN)),
        ),
        _ => format!(
            "UNBOUND ({:?}): e = {:.5}, i = {:.3}°, q = {}, v_inf = {}, ε = {:+.4e} J/kg",
            el.kind,
            e.e,
            e.inclination.to_degrees(),
            u.len(el.periapsis),
            u.speed(el.v_infinity.unwrap_or(0.0)),
            el.specific_energy
        ),
    }
}

/// Osculating elements of every body relative to `primary` (default: the most massive).
pub fn elements_table(system: &NBodySystem, primary: Option<&str>, u: &Units) -> String {
    let Some(p) = (match primary {
        Some(name) => system.find(name),
        None => system.most_massive(),
    }) else {
        return "no bodies".to_string();
    };
    let mut out = String::new();
    let _ = writeln!(out, "Osculating elements relative to {} (two-body μ = G(M + m)):", p.name);
    let _ = writeln!(
        out,
        "{:<14} {:>14} {:>9} {:>8} {:>14} {:>14} {:>14}  bound",
        "name",
        format!("a [{}]", u.spec.length),
        "e",
        "i [deg]",
        format!("T [{}]", u.spec.time),
        format!("q [{}]", u.spec.length),
        format!("Q [{}]", u.spec.length)
    );
    for b in system.bodies() {
        if b.id() == p.id() {
            continue;
        }
        let mu = system.g * (p.mass() + b.mass());
        let el = orbit::elements_from_state(mu, b.position - p.position, b.velocity - p.velocity);
        let _ = writeln!(
            out,
            "{:<14} {:>14} {:>9.5} {:>8.3} {:>14} {:>14} {:>14}  {}",
            truncate(&b.name, 14),
            sig(el.elements.a / u.sys.length),
            el.elements.e,
            el.elements.inclination.to_degrees(),
            el.period.map(|t| sig(t / u.sys.time)).unwrap_or_else(|| "-".into()),
            sig(el.periapsis / u.sys.length),
            el.apoapsis.map(|q| sig(q / u.sys.length)).unwrap_or_else(|| "-".into()),
            if el.is_bound() { "yes" } else { "NO" }
        );
    }
    out
}

pub fn status(sim: &mut Simulation, u: &Units) -> String {
    let d = sim.diagnostics();
    let e_err = sim.conservation.record(&d);
    let p_err = sim.conservation.momentum_error(&d);
    let l_err = sim.conservation.angular_momentum_error(&d);
    let dt = match sim.config().timestep {
        TimestepPolicy::Fixed { dt } => format!("fixed {}", u.time(dt)),
        TimestepPolicy::Adaptive { eta, .. } => format!("adaptive (η={eta}), last {}", u.time(sim.last_dt)),
    };
    let mut out = String::new();
    let _ = writeln!(out, "t = {}   bodies = {}   steps = {}   force evaluations = {}", u.time(d.time), d.bodies, sim.steps, sim.force_evaluations);
    let _ = writeln!(
        out,
        "integrator = {} ({}), dt = {}, forces = {}, collisions = {}",
        sim.config().integrator,
        if sim.config().integrator.symplectic() { "symplectic" } else { "non-symplectic" },
        dt,
        sim.config().force.name(),
        sim.config().collisions
    );
    let _ = writeln!(out, "E = {:+.6e} J  (K = {:.6e}, U = {:+.6e}), virial 2K/|U| = {:.4}", d.total, d.kinetic, d.potential, d.virial_ratio);
    let _ = writeln!(
        out,
        "relative drift: energy {:+.3e} (max {:.3e}), momentum {:.3e}, angular momentum {:.3e}",
        e_err, sim.conservation.max_energy_error, p_err, l_err
    );
    if sim.conservation.collision_energy_change != 0.0 {
        let _ = writeln!(out, "mechanical energy changed by collisions: {:+.4e} J ({} events)", sim.conservation.collision_energy_change, sim.collisions.len());
    }
    out
}

pub fn outcome_label(o: &Outcome, u: &Units) -> String {
    match o {
        Outcome::InFlight { phase } => format!(
            "in flight ({})",
            match phase {
                Phase::Approaching => "approaching",
                Phase::Receding => "receding, unbound",
                Phase::Interacting => "falling back in, unbound",
                Phase::Encounter => "close encounter, verdict suspended",
            }
        ),
        Outcome::Captured { orbits, elements } => format!(
            "CAPTURED into a bound orbit: a = {}, e = {:.4}, period {} ({} full orbit{} completed so far)",
            u.len(elements.elements.a),
            elements.elements.e,
            u.time(elements.period.unwrap_or(f64::NAN)),
            orbits,
            if *orbits == 1 { "" } else { "s" }
        ),
        Outcome::Ejected { v_infinity_out, deflection, specific_energy_change } => format!(
            "EJECTED: leaves at v_inf = {}, deflected by {:.2}°, specific energy change {:+.4e} J/kg",
            u.speed(*v_infinity_out),
            deflection.to_degrees(),
            specific_energy_change
        ),
        Outcome::Merged { into, time } => format!("MERGED into {into} at t = {}", u.time(*time)),
        Outcome::Undetermined => "undetermined".to_string(),
    }
}

pub fn tracker_summary(t: &ProjectileTracker, u: &Units) -> String {
    let mut s = format!("{}: {}", t.name, outcome_label(&t.outcome, u));
    let _ = write!(
        s,
        "\n  closest approach to the barycentre {} at t = {}; periapsis passages {}, apoapsis passages {}",
        u.len(t.closest_approach),
        u.time(t.closest_approach_time),
        t.periapsis_passages,
        t.apoapsis_passages
    );
    if let Some((name, gap)) = &t.closest_body {
        let _ = write!(s, "\n  closest surface-to-surface gap: {} to {name}", u.len(*gap));
    }
    if let Some(name) = &t.encounter_with {
        let _ = write!(s, "\n  currently inside the Hill sphere of {name}");
    }
    let _ = write!(
        s,
        "\n  specific energy relative to the system: {:+.4e} J/kg now, {:+.4e} J/kg at injection",
        t.specific_energy, t.initial_specific_energy
    );
    if !t.absorbed.is_empty() {
        let _ = write!(s, "\n  absorbed: {}", t.absorbed.join(", "));
    }
    if !t.displaced.is_empty() {
        let _ = write!(s, "\n  unbound by the encounter: {}", t.displaced.join(", "));
    }
    s
}

pub fn collisions(sim: &Simulation, u: &Units) -> String {
    if sim.collisions.is_empty() {
        return "no collisions so far\n".to_string();
    }
    let mut out = String::new();
    for e in &sim.collisions {
        let _ = writeln!(out, "t = {:>14}: {}", u.time(e.time), e.describe());
    }
    out
}
