//! `nbody`: command-line driver for the N-body simulator.

mod repl;
mod report;

use anyhow::{anyhow, Context, Result};
use clap::{Parser, Subcommand};
use nbody_core::collision::CollisionMode;
use nbody_core::integrator::IntegratorKind;
use nbody_core::material::Material;
use nbody_core::orbit;
use nbody_core::presets;
use nbody_core::projectile::Outcome;
use nbody_core::scenario::Scenario;
use nbody_core::simulation::{Simulation, TimestepPolicy};
use repl::Shell;
use report::{sig, Units};
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "nbody",
    version,
    about = "Newtonian N-body simulator: symplectic integrators, layered bodies, collisions and projectile capture/ejection analysis",
    long_about = None
)]
struct Cli {
    #[command(subcommand)]
    cmd: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Integrate a scenario file or preset for a duration and report conservation and orbits.
    Run {
        /// Scenario file (.toml) or preset name (see `nbody preset`).
        scenario: String,
        /// Duration in the scenario's time units (default: the scenario's own duration).
        #[arg(long)]
        duration: Option<f64>,
        /// Fixed time step in the scenario's time units.
        #[arg(long)]
        dt: Option<f64>,
        /// leapfrog | yoshida4 | pefrl | rk4
        #[arg(long)]
        integrator: Option<String>,
        /// none | merge | bounce | physical
        #[arg(long)]
        collisions: Option<String>,
        /// Write a trajectory CSV (time, name, x, y, z, vx, vy, vz in scenario units).
        #[arg(long, short)]
        output: Option<PathBuf>,
        /// Number of samples written to the CSV / progress lines.
        #[arg(long, default_value_t = 200)]
        samples: usize,
        /// Save the final state as a scenario file.
        #[arg(long)]
        save: Option<PathBuf>,
    },
    /// Fire a projectile at the system and classify the outcome (captured / ejected / merged).
    Projectile {
        scenario: String,
        /// Excess speed at infinity (length/time units).
        #[arg(long)]
        v_inf: Option<f64>,
        /// Impact parameter (length units).
        #[arg(long)]
        b: Option<f64>,
        /// Projectile mass (mass units).
        #[arg(long)]
        mass: Option<f64>,
        /// Projectile radius (length units).
        #[arg(long)]
        radius: Option<f64>,
        /// Material name (iron, rock, ice, water, gas, plasma, carbon, regolith, rubber).
        #[arg(long)]
        material: Option<String>,
        /// Core material for a two-layer fit (needs --mass, --radius and --material).
        #[arg(long)]
        core: Option<String>,
        /// Approach direction x,y,z (velocity direction at infinity).
        #[arg(long)]
        direction: Option<String>,
        /// Direction of the impact-parameter offset x,y,z.
        #[arg(long)]
        offset: Option<String>,
        /// Start distance from the barycentre (length units).
        #[arg(long)]
        distance: Option<f64>,
        /// Give up after this much simulated time (time units).
        #[arg(long)]
        max_time: Option<f64>,
        #[arg(long)]
        name: Option<String>,
        /// Write a trajectory CSV.
        #[arg(long, short)]
        output: Option<PathBuf>,
        #[arg(long, default_value_t = 400)]
        samples: usize,
    },
    /// Map outcomes over a grid of impact parameters and excess speeds.
    Sweep {
        scenario: String,
        /// Impact parameters as lo:hi:n (length units).
        #[arg(long)]
        b: String,
        /// Excess speeds as lo:hi:n (length/time units).
        #[arg(long)]
        v_inf: String,
        /// Geometric spacing for b.
        #[arg(long)]
        log_b: bool,
        #[arg(long)]
        max_time: Option<f64>,
        #[arg(long)]
        mass: Option<f64>,
        #[arg(long)]
        radius: Option<f64>,
        #[arg(long)]
        material: Option<String>,
        /// Write the grid as JSON.
        #[arg(long)]
        json: Option<PathBuf>,
    },
    /// List presets, or print / save one as a scenario file.
    Preset {
        name: Option<String>,
        #[arg(long, short)]
        output: Option<PathBuf>,
    },
    /// List the built-in materials.
    Materials,
    /// Describe a scenario: bodies, orbital elements, Hill and Roche limits.
    Info { scenario: String },
    /// Interactive shell (add/edit bodies, step, launch projectiles).
    Repl { scenario: Option<String> },
}

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Command::Run { scenario, duration, dt, integrator, collisions, output, samples, save } => {
            cmd_run(&scenario, duration, dt, integrator, collisions, output, samples, save)
        }
        Command::Projectile { scenario, v_inf, b, mass, radius, material, core, direction, offset, distance, max_time, name, output, samples } => {
            let mut kv = HashMap::new();
            let mut put = |k: &str, v: Option<String>| {
                if let Some(v) = v {
                    kv.insert(k.to_string(), v);
                }
            };
            put("vinf", v_inf.map(|x| x.to_string()));
            put("b", b.map(|x| x.to_string()));
            put("mass", mass.map(|x| x.to_string()));
            put("radius", radius.map(|x| x.to_string()));
            put("material", material);
            put("core", core);
            put("dir", direction);
            put("offset", offset);
            put("dist", distance.map(|x| x.to_string()));
            put("name", name);
            cmd_projectile(&scenario, kv, max_time, output, samples)
        }
        Command::Sweep { scenario, b, v_inf, log_b, max_time, mass, radius, material, json } => {
            let mut kv = HashMap::new();
            if let Some(m) = mass {
                kv.insert("mass".into(), m.to_string());
            }
            if let Some(r) = radius {
                kv.insert("radius".into(), r.to_string());
            }
            if let Some(m) = material {
                kv.insert("material".into(), m);
            }
            cmd_sweep(&scenario, &b, &v_inf, log_b, max_time, kv, json)
        }
        Command::Preset { name, output } => cmd_preset(name, output),
        Command::Materials => {
            println!("{:<10} {:>14} {:>12}  description", "name", "density kg/m^3", "restitution");
            for m in Material::BUILTIN.iter() {
                println!("{:<10} {:>14.0} {:>12.2}  {}", m.name(), m.density(), m.restitution(), m.description());
            }
            Ok(())
        }
        Command::Info { scenario } => cmd_info(&scenario),
        Command::Repl { scenario } => {
            let (sc, name) = match scenario {
                Some(s) => {
                    let sc = Shell::load_scenario(&s)?;
                    let name = if sc.name.is_empty() { s } else { sc.name.clone() };
                    (sc, name)
                }
                None => (Scenario { units: nbody_core::scenario::UnitSpec::astronomical(), ..Default::default() }, "empty system".to_string()),
            };
            let mut shell = Shell::from_built(sc.build()?, name)?;
            shell.run()
        }
    }
}

fn open(scenario: &str) -> Result<(Shell, Scenario)> {
    let sc = Shell::load_scenario(scenario)?;
    let name = if sc.name.is_empty() { scenario.to_string() } else { sc.name.clone() };
    let shell = Shell::from_built(sc.build()?, name)?;
    Ok((shell, sc))
}

struct CsvWriter {
    file: std::io::BufWriter<std::fs::File>,
}

impl CsvWriter {
    fn create(path: &PathBuf) -> Result<CsvWriter> {
        let mut file = std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?);
        writeln!(file, "time,name,x,y,z,vx,vy,vz")?;
        Ok(CsvWriter { file })
    }

    fn sample(&mut self, sim: &Simulation, u: &Units) -> Result<()> {
        let t = sim.system.time / u.sys.time;
        for b in sim.system.bodies() {
            let p = b.position / u.sys.length;
            let v = b.velocity / u.sys.velocity();
            writeln!(self.file, "{t:.10e},{},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e},{:.10e}", b.name, p.x, p.y, p.z, v.x, v.y, v.z)?;
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn cmd_run(
    scenario: &str,
    duration: Option<f64>,
    dt: Option<f64>,
    integrator: Option<String>,
    collisions: Option<String>,
    output: Option<PathBuf>,
    samples: usize,
    save: Option<PathBuf>,
) -> Result<()> {
    let (mut shell, sc) = open(scenario)?;
    let u = shell.units.clone();
    if let Some(dt) = dt {
        shell.sim.set_timestep(TimestepPolicy::Fixed { dt: dt * u.sys.time });
    }
    if let Some(i) = integrator {
        shell.sim.set_integrator(IntegratorKind::parse(&i).ok_or_else(|| anyhow!("unknown integrator '{i}'"))?);
    }
    if let Some(c) = collisions {
        shell.sim.set_collisions(CollisionMode::parse(&c).ok_or_else(|| anyhow!("unknown collision mode '{c}'"))?);
    }
    let duration = duration.map(|d| d * u.sys.time).or_else(|| sc.simulation.duration.map(|d| d * u.sys.time)).ok_or_else(|| anyhow!("give --duration (the scenario has none)"))?;

    println!("== {} ==", shell.name);
    if !sc.description.is_empty() {
        println!("{}", sc.description);
    }
    print!("{}", report::body_table(&shell.sim.system, &u));
    print!("{}", report::status(&mut shell.sim, &u));
    println!("running for {} in {} samples", u.time(duration), samples);

    let mut csv = output.as_ref().map(CsvWriter::create).transpose()?;
    if let Some(c) = csv.as_mut() {
        c.sample(&shell.sim, &u)?;
    }
    let t0 = shell.sim.system.time;
    let started = std::time::Instant::now();
    let samples = samples.max(1);
    for k in 1..=samples {
        shell.sim.run_until(t0 + duration * k as f64 / samples as f64);
        if let Some(c) = csv.as_mut() {
            c.sample(&shell.sim, &u)?;
        }
        if k % (samples / 20).max(1) == 0 || k == samples {
            let d = shell.sim.diagnostics();
            let e = shell.sim.conservation.record(&d);
            repl::progress_line(&format!("  t = {:<22} bodies = {:<4} ΔE/E = {:+.2e}   ", u.time(shell.sim.system.time), d.bodies, e));
        }
    }
    repl::progress_end();
    println!("done in {:.2} s wall-clock, {} steps", started.elapsed().as_secs_f64(), shell.sim.steps);
    print!("{}", report::status(&mut shell.sim, &u));
    print!("{}", report::collisions(&shell.sim, &u));
    print!("{}", report::elements_table(&shell.sim.system, None, &u));
    if let Some(path) = save {
        Scenario::snapshot(&shell.sim.system, shell.sim.config(), &u.spec, &shell.name)?.save(&path)?;
        println!("saved final state to {}", path.display());
    }
    Ok(())
}

fn cmd_projectile(scenario: &str, kv: HashMap<String, String>, max_time: Option<f64>, output: Option<PathBuf>, samples: usize) -> Result<()> {
    let (mut shell, sc) = open(scenario)?;
    let u = shell.units.clone();
    let spec = shell.projectile_from_args(&kv)?;
    let max_time = max_time.map(|t| t * u.sys.time).or(shell.projectile_max_time).unwrap_or(100.0 * u.sys.time);

    println!("== {} ==", shell.name);
    if !sc.description.is_empty() {
        println!("{}", sc.description);
    }
    let m_total = shell.sim.system.total_mass();
    let mu = shell.sim.system.g * (m_total + spec.structure.mass());
    println!(
        "projectile {}: mass {}, radius {}, v_inf = {}, b = {}",
        spec.name,
        u.mass(spec.structure.mass()),
        u.len(spec.structure.radius()),
        u.speed(spec.v_infinity),
        u.len(spec.impact_parameter)
    );
    println!(
        "two-body expectation about the barycentre: periapsis {}, deflection {:.2}°",
        u.len(orbit::hyperbolic_periapsis(mu, spec.v_infinity, spec.impact_parameter)),
        orbit::hyperbolic_deflection(mu, spec.v_infinity, spec.impact_parameter.max(1e-300)).to_degrees()
    );
    let id = shell.sim.inject_projectile(&spec)?;
    {
        let b = shell.sim.system.get(id).unwrap();
        println!("injected at {} with {} (start distance {})", u.vec_len(b.position), u.vec_speed(b.velocity), u.len(shell.sim.trackers.last().unwrap().start_distance));
    }
    let mut csv = output.as_ref().map(CsvWriter::create).transpose()?;
    if let Some(c) = csv.as_mut() {
        c.sample(&shell.sim, &u)?;
    }
    let t0 = shell.sim.system.time;
    let samples = samples.max(1);
    let mut last_label = String::new();
    for k in 1..=samples {
        let t_target = t0 + max_time * k as f64 / samples as f64;
        while shell.sim.system.time < t_target {
            let dt = shell.sim.choose_dt().min(t_target - shell.sim.system.time);
            shell.sim.step_with(dt);
            if shell.sim.trackers.last().map(|t| t.outcome.is_terminal()).unwrap_or(true) {
                break;
            }
        }
        if let Some(c) = csv.as_mut() {
            c.sample(&shell.sim, &u)?;
        }
        let tracker = shell.sim.trackers.last().unwrap();
        let label = report::outcome_label(&tracker.outcome, &u);
        if label != last_label {
            println!("  t = {:<22} {}", u.time(shell.sim.system.time), label);
            last_label = label;
        }
        if tracker.outcome.is_terminal() {
            break;
        }
    }
    let tracker = shell.sim.trackers.last().unwrap().clone();
    println!();
    println!("{}", report::tracker_summary(&tracker, &u));
    match &tracker.outcome {
        Outcome::Captured { elements, .. } => println!("  orbit: {}", report::elements_line(elements, &u)),
        Outcome::InFlight { .. } | Outcome::Undetermined => {
            println!("  no verdict within {} of simulated time; increase --max-time", u.time(max_time));
            if let Some(el) = tracker.current_elements(&shell.sim.system) {
                println!("  current osculating orbit: {}", report::elements_line(&el, &u));
            }
        }
        _ => {}
    }
    print!("{}", report::collisions(&shell.sim, &u));
    print!("{}", report::status(&mut shell.sim, &u));
    Ok(())
}

fn cmd_sweep(scenario: &str, b: &str, v_inf: &str, log_b: bool, max_time: Option<f64>, kv: HashMap<String, String>, json: Option<PathBuf>) -> Result<()> {
    let (shell, _) = open(scenario)?;
    let u = shell.units.clone();
    let template = shell.projectile_from_args(&kv).or_else(|e| shell.projectile.clone().ok_or(e))?;
    let bs: Vec<f64> = repl::parse_range(b, log_b)?.into_iter().map(|x| x * u.sys.length).collect();
    let vs: Vec<f64> = repl::parse_range(v_inf, false)?.into_iter().map(|x| x * u.sys.velocity()).collect();
    let max_time = max_time.map(|t| t * u.sys.time).or(shell.projectile_max_time).unwrap_or(100.0 * u.sys.time);
    println!("== {} == sweep of {} × {} projectile experiments, max {} each", shell.name, vs.len(), bs.len(), u.time(max_time));
    let result = shell.run_sweep(&template, &bs, &vs, max_time);
    print!("{}", shell.sweep_report(&result));
    if let Some(path) = json {
        std::fs::write(&path, serde_json::to_string_pretty(&result)?)?;
        println!("wrote {}", path.display());
    }
    Ok(())
}

fn cmd_preset(name: Option<String>, output: Option<PathBuf>) -> Result<()> {
    match name {
        None => {
            for (n, d) in presets::names() {
                println!("{n:<20} {d}");
            }
            Ok(())
        }
        Some(n) => {
            let sc = presets::preset(&n).ok_or_else(|| anyhow!("unknown preset '{n}'"))?;
            match output {
                Some(path) => {
                    sc.save(&path)?;
                    println!("wrote {}", path.display());
                }
                None => print!("{}", sc.to_toml()?),
            }
            Ok(())
        }
    }
}

fn cmd_info(scenario: &str) -> Result<()> {
    let (mut shell, sc) = open(scenario)?;
    let u = shell.units.clone();
    println!("== {} ==", shell.name);
    if !sc.description.is_empty() {
        println!("{}", sc.description);
    }
    println!("units: {}", u.sys);
    print!("{}", report::body_table(&shell.sim.system, &u));
    print!("{}", report::status(&mut shell.sim, &u));
    print!("{}", report::elements_table(&shell.sim.system, None, &u));
    if let Some(primary) = shell.sim.system.most_massive().cloned() {
        println!("Stability scales relative to {}:", primary.name);
        println!("{:<14} {:>16} {:>16} {:>16}  {}", "name", format!("Hill r [{}]", u.spec.length), "Roche rigid", "Roche fluid", "note");
        for b in shell.sim.system.bodies() {
            if b.id() == primary.id() {
                continue;
            }
            let mu = shell.sim.system.g * (primary.mass() + b.mass());
            let el = orbit::elements_from_state(mu, b.position - primary.position, b.velocity - primary.velocity);
            let hill = orbit::hill_radius(el.elements.a.abs(), el.elements.e.min(0.999), b.mass(), primary.mass());
            let rr = orbit::roche_limit_rigid(primary.radius(), primary.mean_density(), b.mean_density());
            let rf = orbit::roche_limit_fluid(primary.radius(), primary.mean_density(), b.mean_density());
            let d = (b.position - primary.position).norm();
            let note = if d < rr {
                "inside rigid Roche limit!"
            } else if d < rf {
                "inside fluid Roche limit"
            } else if !el.is_bound() {
                "unbound"
            } else {
                ""
            };
            println!(
                "{:<14} {:>16} {:>16} {:>16}  {}",
                b.name,
                sig(hill / u.sys.length),
                sig(rr / u.sys.length),
                sig(rf / u.sys.length),
                note
            );
        }
    }
    let t_dyn = nbody_core::force::min_dynamical_time(&shell.sim.system.positions(), &shell.sim.system.masses(), shell.sim.system.g, shell.sim.system.softening);
    if t_dyn.is_finite() {
        println!(
            "shortest two-body dynamical time {} (period {}); a fixed step of {} gives 1000 steps per shortest orbit",
            u.time(t_dyn),
            u.time(2.0 * std::f64::consts::PI * t_dyn),
            u.time(2.0 * std::f64::consts::PI * t_dyn / 1000.0)
        );
    }
    if let Some(p) = &shell.projectile {
        println!(
            "projectile: {} ({}), v_inf = {}, b = {}, direction {}",
            p.name,
            u.mass(p.structure.mass()),
            u.speed(p.v_infinity),
            u.len(p.impact_parameter),
            format_args!("[{}, {}, {}]", sig(p.direction.x), sig(p.direction.y), sig(p.direction.z))
        );
    }
    Ok(())
}
