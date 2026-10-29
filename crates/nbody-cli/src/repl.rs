//! Interactive shell: edit the system, step it, launch projectiles.

use crate::report::{self, sig, Units};
use anyhow::{anyhow, bail, Context, Result};
use nbody_core::body::{Body, BodyKind};
use nbody_core::collision::CollisionMode;
use nbody_core::force::ForceMethod;
use nbody_core::integrator::IntegratorKind;
use nbody_core::material::Material;
use nbody_core::orbit::Elements;
use nbody_core::presets;
use nbody_core::projectile::{self, ProjectileSpec};
use nbody_core::scenario::{Built, Scenario};
use nbody_core::simulation::{Simulation, TimestepPolicy};
use nbody_core::structure::{self, SphericalStructure};
use nbody_core::vector::Vec3;
use std::collections::HashMap;
use std::io::{BufRead, Write};

pub const HELP: &str = r#"Commands (quantities in the scenario's units; angles in degrees):
  list | ls                         table of bodies
  show <name>                       structure, orbit, Hill and Roche limits of one body
  elements [primary]                osculating elements of every body
  status                            time, energy drift, integrator settings
  collisions                        collision log

  add <name> mass=<m> [material=<mat>] [radius=<r>] [core=<mat>] [kind=<k>]
             [pos=x,y,z] [vel=x,y,z] | [orbit=<around>,a[,e[,i[,raan[,argp[,nu]]]]]]
  remove <name>
  set <name> mass=<m> | radius=<r> | material=<mat> | pos=x,y,z | vel=x,y,z | name=<n> | kind=<k>
  shell <name> add <thickness> <material>     add an outer shell
  shell <name> pop                            remove the outer shell
  circular <name> around <other>              circular orbit at the current distance
  orbit <name> around <other> a=<a> [e=] [i=] [raan=] [argp=] [nu=]
  barycentre                                  move to the centre-of-mass frame

  step [n]                          advance n steps (default 1)
  run <duration>                    advance by a duration
  dt <value> | adaptive <eta>[,<dt_min>,<dt_max>]
  integrator leapfrog|yoshida4|pefrl|rk4
  collisions none|merge|bounce|physical
  softening <length> | force direct|bh[,theta]

  projectile vinf=<v> b=<b> (mass=<m>|radius=<r>) [material=] [core=] [name=] [dir=x,y,z] [offset=x,y,z] [dist=<d>]
  outcome                           status of the latest projectile
  sweep vinf=<lo>:<hi>:<n> b=<lo>:<hi>:<n> [maxtime=<t>] [log] [projectile options]

  materials | presets
  save <file.toml> | load <file.toml | preset-name>
  units <length> <mass> <time>      change display/input units (e.g. units au solar year)
  help | quit
"#;

/// Overwrite the current stderr line on a terminal; print a new line otherwise.
pub fn progress_line(text: &str) {
    use std::io::IsTerminal;
    if std::io::stderr().is_terminal() {
        eprint!("\r{text}");
    } else {
        eprintln!("{}", text.trim_end());
    }
}

pub fn progress_end() {
    use std::io::IsTerminal;
    if std::io::stderr().is_terminal() {
        eprintln!();
    }
}

pub struct Shell {
    pub sim: Simulation,
    pub units: Units,
    pub name: String,
    pub projectile: Option<ProjectileSpec>,
    pub projectile_max_time: Option<f64>,
    pub quiet: bool,
}

fn parse_kv(tokens: &[&str]) -> (Vec<String>, HashMap<String, String>) {
    let mut positional = Vec::new();
    let mut kv = HashMap::new();
    for t in tokens {
        if let Some((k, v)) = t.split_once('=') {
            kv.insert(k.trim().to_ascii_lowercase(), v.trim().to_string());
        } else {
            positional.push(t.to_string());
        }
    }
    (positional, kv)
}

fn parse_f64(s: &str, what: &str) -> Result<f64> {
    s.parse::<f64>().with_context(|| format!("{what}: '{s}' is not a number"))
}

fn parse_vec3(s: &str, what: &str) -> Result<Vec3> {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() != 3 {
        bail!("{what}: expected x,y,z, got '{s}'");
    }
    Ok(Vec3::new(parse_f64(parts[0], what)?, parse_f64(parts[1], what)?, parse_f64(parts[2], what)?))
}

fn parse_material(s: &str) -> Result<Material> {
    Material::parse(s).ok_or_else(|| {
        anyhow!("unknown material '{s}'; known: {}", Material::BUILTIN.iter().map(|m| m.name()).collect::<Vec<_>>().join(", "))
    })
}

/// `lo:hi:n` → n values, linear or geometric.
pub fn parse_range(s: &str, log: bool) -> Result<Vec<f64>> {
    let parts: Vec<&str> = s.split(':').collect();
    match parts.len() {
        1 => Ok(vec![parse_f64(parts[0], "range")?]),
        3 => {
            let lo = parse_f64(parts[0], "range start")?;
            let hi = parse_f64(parts[1], "range end")?;
            let n: usize = parts[2].parse().context("range count")?;
            if n == 0 {
                bail!("range count must be positive");
            }
            if n == 1 {
                return Ok(vec![lo]);
            }
            Ok((0..n)
                .map(|k| {
                    let t = k as f64 / (n - 1) as f64;
                    if log && lo > 0.0 && hi > 0.0 {
                        lo * (hi / lo).powf(t)
                    } else {
                        lo + (hi - lo) * t
                    }
                })
                .collect())
        }
        _ => bail!("expected lo:hi:n, got '{s}'"),
    }
}

impl Shell {
    pub fn from_built(built: Built, name: String) -> Result<Shell> {
        let units = Units::new(built.unit_spec.clone())?;
        Ok(Shell {
            sim: Simulation::new(built.system, built.config),
            units,
            name,
            projectile: built.projectile,
            projectile_max_time: built.projectile_max_time,
            quiet: false,
        })
    }

    pub fn load_scenario(arg: &str) -> Result<Scenario> {
        let path = std::path::Path::new(arg);
        if path.exists() {
            return Ok(Scenario::load(path)?);
        }
        if let Some(p) = presets::preset(arg) {
            return Ok(p);
        }
        bail!(
            "'{arg}' is neither a file nor a preset; presets: {}",
            presets::names().iter().map(|(n, _)| *n).collect::<Vec<_>>().join(", ")
        )
    }

    /// Build a structure from REPL-style arguments.
    pub fn structure_from_args(&self, kv: &HashMap<String, String>) -> Result<SphericalStructure> {
        let u = &self.units.sys;
        let mass = kv.get("mass").map(|s| parse_f64(s, "mass")).transpose()?.map(|m| m * u.mass);
        let radius = kv.get("radius").map(|s| parse_f64(s, "radius")).transpose()?.map(|r| r * u.length);
        let material = kv.get("material").map(|s| parse_material(s)).transpose()?;
        let core = kv.get("core").map(|s| parse_material(s)).transpose()?;
        Ok(match (mass, radius) {
            (Some(m), None) => SphericalStructure::uniform_with_mass(material.unwrap_or(Material::Rock), m)?,
            (None, Some(r)) => SphericalStructure::uniform_with_radius(material.unwrap_or(Material::Rock), r)?,
            (Some(m), Some(r)) => {
                if let (Some(c), Some(mantle)) = (core, material.clone()) {
                    structure::fit_two_layer(c, mantle, m, r)?
                } else {
                    let base = material.unwrap_or(Material::Rock);
                    let density = m / structure::sphere_volume(r);
                    SphericalStructure::uniform_with_radius(Material::custom(base.name(), density, base.restitution()), r)?
                }
            }
            (None, None) => bail!("give mass=<m> and/or radius=<r>"),
        })
    }

    pub fn run(&mut self) -> Result<()> {
        println!("N-body shell: '{}' with {} bodies. Type 'help' for commands.", self.name, self.sim.system.len());
        print!("{}", report::body_table(&self.sim.system, &self.units));
        let stdin = std::io::stdin();
        let mut out = std::io::stdout();
        loop {
            print!("nbody> ");
            out.flush()?;
            let mut line = String::new();
            if stdin.lock().read_line(&mut line)? == 0 {
                println!();
                break;
            }
            match self.execute(line.trim()) {
                Ok(true) => {}
                Ok(false) => break,
                Err(e) => println!("error: {e:#}"),
            }
        }
        Ok(())
    }

    fn body_id(&self, name: &str) -> Result<nbody_core::body::BodyId> {
        Ok(self.sim.system.id_of(name)?)
    }

    pub fn execute(&mut self, line: &str) -> Result<bool> {
        if line.is_empty() || line.starts_with('#') {
            return Ok(true);
        }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        let cmd = tokens[0].to_ascii_lowercase();
        let args = &tokens[1..];
        let u = self.units.sys;
        match cmd.as_str() {
            "help" | "?" => print!("{HELP}"),
            "quit" | "exit" | "q" => return Ok(false),
            "list" | "ls" => print!("{}", report::body_table(&self.sim.system, &self.units)),
            "show" => {
                let name = args.first().ok_or_else(|| anyhow!("show <name>"))?;
                let b = self.sim.system.find(name).ok_or_else(|| anyhow!("no body '{name}'"))?;
                print!("{}", report::body_detail(&self.sim.system, b, &self.units));
            }
            "elements" => print!("{}", report::elements_table(&self.sim.system, args.first().copied(), &self.units)),
            "status" | "energy" => print!("{}", report::status(&mut self.sim, &self.units)),
            "collisions" if args.is_empty() => print!("{}", report::collisions(&self.sim, &self.units)),
            "materials" => {
                for m in Material::BUILTIN.iter() {
                    println!("  {:<9} rho = {:>5.0} kg/m^3  e = {:.2}  {}", m.name(), m.density(), m.restitution(), m.description());
                }
            }
            "presets" => {
                for (n, d) in presets::names() {
                    println!("  {n:<20} {d}");
                }
            }
            "add" => {
                let (pos, kv) = parse_kv(args);
                let name = pos.first().ok_or_else(|| anyhow!("add <name> mass=<m> ..."))?.clone();
                if self.sim.system.find(&name).is_some() {
                    bail!("a body named '{name}' already exists");
                }
                let structure = self.structure_from_args(&kv)?;
                let kind = kv.get("kind").map(|k| BodyKind::parse(k).ok_or_else(|| anyhow!("unknown kind '{k}'"))).transpose()?.unwrap_or_default();
                let position = kv.get("pos").map(|s| parse_vec3(s, "pos")).transpose()?.unwrap_or(Vec3::ZERO) * u.length;
                let velocity = kv.get("vel").map(|s| parse_vec3(s, "vel")).transpose()?.unwrap_or(Vec3::ZERO) * u.velocity();
                let id = self.sim.system.add(Body::new(name.clone(), structure, position, velocity).with_kind(kind));
                if let Some(o) = kv.get("orbit") {
                    let parts: Vec<&str> = o.split(',').collect();
                    let around = self.body_id(parts[0]).map_err(|_| anyhow!("orbit: no body '{}'", parts[0]))?;
                    let num = |i: usize, default: f64| -> Result<f64> { parts.get(i).map(|s| parse_f64(s, "orbit")).transpose().map(|v| v.unwrap_or(default)) };
                    let a = num(1, 1.0)? * u.length;
                    let el = Elements {
                        a,
                        e: num(2, 0.0)?,
                        inclination: num(3, 0.0)?.to_radians(),
                        raan: num(4, 0.0)?.to_radians(),
                        arg_periapsis: num(5, 0.0)?.to_radians(),
                        true_anomaly: num(6, 0.0)?.to_radians(),
                    };
                    self.sim.system.set_orbit(id, around, &el)?;
                }
                self.sim.reset_baseline();
                let b = self.sim.system.get(id).unwrap();
                println!("added {} ({}, R = {})", b.name, self.units.mass(b.mass()), self.units.len(b.radius()));
            }
            "remove" | "rm" => {
                let name = args.first().ok_or_else(|| anyhow!("remove <name>"))?;
                let id = self.body_id(name)?;
                self.sim.system.remove(id);
                self.sim.reset_baseline();
                println!("removed {name}");
            }
            "set" => {
                let (pos, kv) = parse_kv(args);
                let name = pos.first().ok_or_else(|| anyhow!("set <name> key=value ..."))?;
                let id = self.body_id(name)?;
                let b = self.sim.system.get_mut(id).unwrap();
                for (k, v) in &kv {
                    match k.as_str() {
                        "mass" => b.set_mass(parse_f64(v, "mass")? * u.mass)?,
                        "radius" => b.set_radius(parse_f64(v, "radius")? * u.length)?,
                        "material" => b.set_material(parse_material(v)?)?,
                        "pos" => b.position = parse_vec3(v, "pos")? * u.length,
                        "vel" => b.velocity = parse_vec3(v, "vel")? * u.velocity(),
                        "name" => b.name = v.clone(),
                        "kind" => b.kind = BodyKind::parse(v).ok_or_else(|| anyhow!("unknown kind '{v}'"))?,
                        other => bail!("unknown property '{other}'"),
                    }
                }
                self.sim.reset_baseline();
                let b = self.sim.system.get(id).unwrap();
                print!("{}", report::body_detail(&self.sim.system, b, &self.units));
            }
            "shell" => {
                let name = args.first().ok_or_else(|| anyhow!("shell <name> add <thickness> <material> | pop"))?;
                let id = self.body_id(name)?;
                let b = self.sim.system.get_mut(id).unwrap();
                let mut s = b.structure().clone();
                match args.get(1).map(|s| s.to_ascii_lowercase()).as_deref() {
                    Some("add") => {
                        let t = parse_f64(args.get(2).ok_or_else(|| anyhow!("thickness missing"))?, "thickness")? * u.length;
                        let m = parse_material(args.get(3).ok_or_else(|| anyhow!("material missing"))?)?;
                        s.push_shell(t, m)?;
                    }
                    Some("pop") => {
                        s.pop_shell().ok_or_else(|| anyhow!("a body must keep at least one layer"))?;
                    }
                    _ => bail!("shell <name> add <thickness> <material> | shell <name> pop"),
                }
                b.set_structure(s);
                self.sim.reset_baseline();
                let b = self.sim.system.get(id).unwrap();
                print!("{}", report::body_detail(&self.sim.system, b, &self.units));
            }
            "circular" => {
                if args.len() < 3 || !args[1].eq_ignore_ascii_case("around") {
                    bail!("circular <name> around <other>");
                }
                let id = self.body_id(args[0])?;
                let around = self.body_id(args[2])?;
                self.sim.system.set_circular_orbit(id, around, Vec3::Z)?;
                self.sim.reset_baseline();
                println!("{}", report::elements_table(&self.sim.system, Some(args[2]), &self.units));
            }
            "orbit" => {
                if args.len() < 3 || !args[1].eq_ignore_ascii_case("around") {
                    bail!("orbit <name> around <other> a=<a> [e= i= raan= argp= nu=]");
                }
                let id = self.body_id(args[0])?;
                let around = self.body_id(args[2])?;
                let (_, kv) = parse_kv(&args[3..]);
                let get = |k: &str, d: f64| -> Result<f64> { kv.get(k).map(|s| parse_f64(s, k)).transpose().map(|v| v.unwrap_or(d)) };
                let el = Elements {
                    a: get("a", 1.0)? * u.length,
                    e: get("e", 0.0)?,
                    inclination: get("i", 0.0)?.to_radians(),
                    raan: get("raan", 0.0)?.to_radians(),
                    arg_periapsis: get("argp", 0.0)?.to_radians(),
                    true_anomaly: get("nu", 0.0)?.to_radians(),
                };
                self.sim.system.set_orbit(id, around, &el)?;
                self.sim.reset_baseline();
                println!("{}", report::elements_table(&self.sim.system, Some(args[2]), &self.units));
            }
            "barycentre" | "barycenter" => {
                self.sim.system.move_to_barycentre();
                self.sim.reset_baseline();
                println!("moved to the centre-of-mass frame");
            }
            "step" => {
                let n: usize = args.first().map(|s| s.parse()).transpose().context("step count")?.unwrap_or(1);
                for _ in 0..n {
                    self.sim.step();
                }
                print!("{}", report::status(&mut self.sim, &self.units));
            }
            "run" => {
                let d = parse_f64(args.first().ok_or_else(|| anyhow!("run <duration>"))?, "duration")? * u.time;
                self.advance_with_progress(d);
                print!("{}", report::status(&mut self.sim, &self.units));
            }
            "dt" => {
                let dt = parse_f64(args.first().ok_or_else(|| anyhow!("dt <value>"))?, "dt")? * u.time;
                self.sim.set_timestep(TimestepPolicy::Fixed { dt });
                println!("fixed step {}", self.units.time(dt));
            }
            "adaptive" => {
                let parts: Vec<&str> = args.first().ok_or_else(|| anyhow!("adaptive <eta>[,<dt_min>,<dt_max>]"))?.split(',').collect();
                let eta = parse_f64(parts[0], "eta")?;
                let dt_min = parts.get(1).map(|s| parse_f64(s, "dt_min")).transpose()?.unwrap_or(0.0) * u.time;
                let dt_max = parts.get(2).map(|s| parse_f64(s, "dt_max")).transpose()?.unwrap_or(f64::INFINITY) * u.time;
                self.sim.set_timestep(TimestepPolicy::Adaptive { eta, dt_min, dt_max });
                println!("adaptive steps: eta = {eta}");
            }
            "integrator" => {
                let k = IntegratorKind::parse(args.first().ok_or_else(|| anyhow!("integrator <kind>"))?).ok_or_else(|| anyhow!("unknown integrator"))?;
                self.sim.set_integrator(k);
                println!("{}: {}", k, k.description());
            }
            "collisions" => {
                let m = CollisionMode::parse(args[0]).ok_or_else(|| anyhow!("unknown collision mode"))?;
                self.sim.set_collisions(m);
                println!("collisions: {m}");
            }
            "softening" => {
                let e = parse_f64(args.first().ok_or_else(|| anyhow!("softening <length>"))?, "softening")? * u.length;
                self.sim.system.softening = e;
                self.sim.reset_baseline();
                println!("softening {}", self.units.len(e));
            }
            "force" => {
                let spec = args.first().ok_or_else(|| anyhow!("force direct|bh[,theta]"))?;
                let parts: Vec<&str> = spec.split(',').collect();
                let f = match parts[0].to_ascii_lowercase().as_str() {
                    "direct" => ForceMethod::Direct,
                    "bh" | "barnes_hut" | "tree" => ForceMethod::BarnesHut { theta: parts.get(1).map(|s| parse_f64(s, "theta")).transpose()?.unwrap_or(0.5) },
                    _ => bail!("force direct|bh[,theta]"),
                };
                self.sim.set_force(f);
                println!("forces: {f:?}");
            }
            "projectile" => {
                let (_, kv) = parse_kv(args);
                let spec = self.projectile_from_args(&kv)?;
                let id = self.sim.inject_projectile(&spec)?;
                let b = self.sim.system.get(id).unwrap();
                println!(
                    "launched {} ({}) from {} at {}; v_inf = {}, b = {}",
                    b.name,
                    self.units.mass(b.mass()),
                    self.units.vec_len(b.position),
                    self.units.vec_speed(b.velocity),
                    self.units.speed(spec.v_infinity),
                    self.units.len(spec.impact_parameter)
                );
                self.projectile = Some(spec);
            }
            "outcome" => {
                let t = self.sim.trackers.last().ok_or_else(|| anyhow!("no projectile launched yet"))?;
                println!("{}", report::tracker_summary(t, &self.units));
                if let Some(el) = t.current_elements(&self.sim.system) {
                    println!("  now: {}", report::elements_line(&el, &self.units));
                }
            }
            "sweep" => {
                let (flags, kv) = parse_kv(args);
                let log = flags.iter().any(|f| f.eq_ignore_ascii_case("log"));
                let template = self.projectile_from_args(&kv).or_else(|e| self.projectile.clone().ok_or(e))?;
                let bs = parse_range(kv.get("b").ok_or_else(|| anyhow!("b=lo:hi:n"))?, log)?.into_iter().map(|b| b * u.length).collect::<Vec<_>>();
                let vs = parse_range(kv.get("vinf").ok_or_else(|| anyhow!("vinf=lo:hi:n"))?, false)?.into_iter().map(|v| v * u.velocity()).collect::<Vec<_>>();
                let max_time = kv.get("maxtime").map(|s| parse_f64(s, "maxtime")).transpose()?.map(|t| t * u.time).or(self.projectile_max_time).unwrap_or(100.0 * u.time);
                let result = self.run_sweep(&template, &bs, &vs, max_time);
                print!("{}", self.sweep_report(&result));
            }
            "save" => {
                let path = args.first().ok_or_else(|| anyhow!("save <file.toml>"))?;
                let sc = Scenario::snapshot(&self.sim.system, self.sim.config(), &self.units.spec, &self.name)?;
                sc.save(path)?;
                println!("saved {} bodies to {path}", self.sim.system.len());
            }
            "load" => {
                let arg = args.first().ok_or_else(|| anyhow!("load <file.toml | preset>"))?;
                let sc = Self::load_scenario(arg)?;
                let name = if sc.name.is_empty() { arg.to_string() } else { sc.name.clone() };
                let built = sc.build()?;
                *self = Shell::from_built(built, name)?;
                print!("{}", report::body_table(&self.sim.system, &self.units));
            }
            "units" => {
                if args.len() != 3 {
                    bail!("units <length> <mass> <time>");
                }
                self.units = Units::new(nbody_core::scenario::UnitSpec::new(args[0], args[1], args[2]))?;
                println!("display units: {}", self.units.sys);
            }
            "time" => println!("t = {}", self.units.time(self.sim.system.time)),
            other => bail!("unknown command '{other}' (try 'help')"),
        }
        Ok(true)
    }

    pub fn projectile_from_args(&self, kv: &HashMap<String, String>) -> Result<ProjectileSpec> {
        let u = self.units.sys;
        let base = self.projectile.clone();
        let structure = if kv.contains_key("mass") || kv.contains_key("radius") {
            self.structure_from_args(kv)?
        } else if let Some(b) = &base {
            b.structure.clone()
        } else {
            bail!("give mass=<m> or radius=<r> for the projectile");
        };
        let v_inf = match kv.get("vinf").or(kv.get("v_inf")) {
            Some(s) => parse_f64(s, "vinf")? * u.velocity(),
            None => base.as_ref().map(|b| b.v_infinity).ok_or_else(|| anyhow!("give vinf=<speed>"))?,
        };
        let b = match kv.get("b") {
            Some(s) => parse_f64(s, "b")? * u.length,
            None => base.as_ref().map(|b| b.impact_parameter).unwrap_or(0.0),
        };
        let mut spec = ProjectileSpec::new(
            kv.get("name").cloned().or_else(|| base.as_ref().map(|b| b.name.clone())).unwrap_or_else(|| "projectile".into()),
            structure,
            v_inf,
            b,
        );
        spec.direction = match kv.get("dir") {
            Some(s) => parse_vec3(s, "dir")?,
            None => base.as_ref().map(|b| b.direction).unwrap_or(Vec3::X),
        };
        spec.offset_direction = match kv.get("offset") {
            Some(s) => Some(parse_vec3(s, "offset")?),
            None => base.as_ref().and_then(|b| b.offset_direction),
        };
        spec.start_distance = match kv.get("dist") {
            Some(s) => Some(parse_f64(s, "dist")? * u.length),
            None => base.as_ref().and_then(|b| b.start_distance),
        };
        Ok(spec)
    }

    pub fn advance_with_progress(&mut self, duration: f64) {
        let t0 = self.sim.system.time;
        let chunks = 20;
        for k in 1..=chunks {
            self.sim.run_until(t0 + duration * k as f64 / chunks as f64);
            if !self.quiet && duration * self.units.sys.time > 0.0 {
                let d = self.sim.diagnostics();
                let e = self.sim.conservation.record(&d);
                progress_line(&format!("  t = {:<22} bodies = {:<4} ΔE/E = {:+.2e}   ", self.units.time(self.sim.system.time), d.bodies, e));
            }
        }
        if !self.quiet {
            progress_end();
        }
    }

    pub fn run_sweep(&self, template: &ProjectileSpec, bs: &[f64], vs: &[f64], max_time: f64) -> projectile::SweepResult {
        let total = bs.len() * vs.len();
        let mut progress = |done: usize, _total: usize| {
            if done % (total / 20).max(1) == 0 || done == total {
                progress_line(&format!("  sweep {done}/{total}   "));
            }
        };
        let r = projectile::sweep(&self.sim.system, self.sim.config(), template, bs, vs, max_time, Some(&mut progress));
        progress_end();
        r
    }

    pub fn sweep_report(&self, r: &projectile::SweepResult) -> String {
        let u = &self.units;
        let mut out = String::new();
        out.push_str(&format!(
            "rows: v_inf [{}/{}], columns: impact parameter b from {} to {} ({} values)\n",
            u.spec.length,
            u.spec.time,
            sig(r.impact_parameters[0] / u.sys.length),
            sig(r.impact_parameters[r.impact_parameters.len() - 1] / u.sys.length),
            r.impact_parameters.len()
        ));
        out.push_str("C = captured   . = ejected   # = merged   ? = undetermined within max time\n");
        let fr = r.cross_section_fractions();
        for (iv, row) in r.grid.iter().enumerate() {
            let line: String = row.iter().map(|c| c.symbol()).collect();
            out.push_str(&format!(
                "{:>12} | {}  | σ-weighted: capture {:>5.1}%  merge {:>5.1}%  eject {:>5.1}%\n",
                sig(r.v_infinities[iv] / u.sys.velocity()),
                line,
                100.0 * fr[iv][0],
                100.0 * fr[iv][2],
                100.0 * fr[iv][1]
            ));
        }
        out
    }
}
