//! `nbody-gui`: an interactive front-end for the N-body simulator.
//!
//! The same [`App`] runs in a desktop window (`run_native`, used by the
//! `nbody-gui` binary) and, compiled to `wasm32-unknown-unknown`, inside a
//! web page (`start`, exported through wasm-bindgen; see `scripts/build.mjs`).
//!
//! Left panel: bodies (select, edit mass / radius / material / layers,
//! place on orbits, add, delete). Right panel: projectile mode and outcome
//! sweeps. Bottom: conservation diagnostics. Centre: the orbit view
//! (drag to pan, scroll to zoom, right-drag to tilt, click to select).

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};
use nbody_core::body::{Body, BodyId, BodyKind};
use nbody_core::collision::CollisionMode;
use nbody_core::integrator::IntegratorKind;
use nbody_core::material::Material;
use nbody_core::orbit::{self, Elements};
use nbody_core::presets;
use nbody_core::projectile::{self, Outcome, Phase, ProjectileSpec};
use nbody_core::scenario::{Scenario, UnitSpec};
use nbody_core::simulation::{Simulation, TimestepPolicy};
use nbody_core::structure::{self, SphericalStructure};
use nbody_core::units::UnitSystem;
use nbody_core::vector::Vec3;
use std::collections::{HashMap, VecDeque};

/// Run as a native desktop window.
#[cfg(not(target_arch = "wasm32"))]
pub fn run_native() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1500.0, 940.0]).with_title("N-body problem"),
        ..Default::default()
    };
    eframe::run_native(
        "N-body problem",
        options,
        Box::new(|cc| {
            cc.egui_ctx.set_visuals(egui::Visuals::dark());
            Ok(Box::new(App::new()))
        }),
    )
}

/// A drag value whose step scales with the magnitude of the value.
fn drag(ui: &mut egui::Ui, v: &mut f64, min_speed: f64) -> egui::Response {
    let speed = v.abs() * 0.02 + min_speed;
    ui.add(egui::DragValue::new(v).speed(speed))
}

fn c32(c: [u8; 3]) -> Color32 {
    Color32::from_rgb(c[0], c[1], c[2])
}

/// Compact number formatting for labels.
fn sig(x: f64) -> String {
    if x == 0.0 {
        return "0".into();
    }
    if !x.is_finite() {
        return format!("{x}");
    }
    let a = x.abs();
    if (1e-3..1e6).contains(&a) {
        let decimals = (4 - a.log10().floor() as i32).max(0) as usize;
        let s = format!("{x:.decimals$}");
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    } else {
        format!("{x:.3e}")
    }
}

#[derive(Clone, Copy)]
struct View {
    center: Vec3,
    /// Pixels per metre.
    scale: f64,
    yaw: f64,
    pitch: f64,
}

impl View {
    fn rotate(&self, p: Vec3) -> Vec3 {
        let d = p - self.center;
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        let x1 = cy * d.x - sy * d.y;
        let y1 = sy * d.x + cy * d.y;
        let z1 = d.z;
        Vec3::new(x1, cp * y1 - sp * z1, sp * y1 + cp * z1)
    }

    fn unrotate(&self, v: Vec3) -> Vec3 {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        let y1 = cp * v.y + sp * v.z;
        let z1 = -sp * v.y + cp * v.z;
        let x1 = v.x;
        Vec3::new(cy * x1 + sy * y1, -sy * x1 + cy * y1, z1)
    }

    fn project(&self, p: Vec3, rect: Rect) -> (Pos2, f64) {
        let v = self.rotate(p);
        (Pos2::new(rect.center().x + (v.x * self.scale) as f32, rect.center().y - (v.y * self.scale) as f32), v.z)
    }

    fn pan(&mut self, delta: Vec2) {
        let d = Vec3::new(-delta.x as f64 / self.scale, delta.y as f64 / self.scale, 0.0);
        self.center += self.unrotate(d);
    }
}

pub struct App {
    sim: Simulation,
    scenario: Scenario,
    scenario_name: String,
    units: UnitSystem,
    unit_spec: UnitSpec,
    projectile_template: Option<ProjectileSpec>,

    running: bool,
    /// Display time units of simulated time per real second.
    time_scale: f64,
    max_steps_per_frame: usize,
    steps_last_frame: usize,
    lagging: bool,

    view: View,
    follow: Option<BodyId>,
    show_labels: bool,
    show_trails: bool,
    trail_len: usize,
    trails: HashMap<BodyId, VecDeque<Vec3>>,

    selected: Option<BodyId>,
    edit_loaded_for: Option<BodyId>,
    edit_name: String,
    shell_thickness: f64,
    shell_material: Material,
    orbit_target: Option<BodyId>,

    new_name: String,
    new_mass: f64,
    new_radius: f64,
    new_material: Material,
    new_core: Option<Material>,
    new_kind: BodyKind,
    new_around: Option<BodyId>,
    new_a: f64,
    new_e: f64,
    new_i: f64,
    new_pos: [f64; 3],
    new_vel: [f64; 3],

    p_name: String,
    p_v_inf: f64,
    p_b: f64,
    p_mass: f64,
    p_material: Material,
    p_angle: f64,
    p_side: bool,
    p_start: f64,
    p_max_time: f64,
    sweep_b: [f64; 2],
    sweep_n: usize,
    sweep_v: [f64; 2],
    sweep_nv: usize,
    sweep_result: String,

    fixed_dt: f64,
    adaptive: bool,
    eta: f64,

    preset_index: usize,
    show_bodies: bool,
    show_projectile: bool,
    layout_initialised: bool,
    #[cfg_attr(target_arch = "wasm32", allow(dead_code))]
    file_path: String,
    energy_history: VecDeque<f32>,
    log: VecDeque<String>,
    error: Option<String>,
}

impl App {
    pub fn new() -> App {
        let scenario = presets::preset("solar-system").expect("preset");
        let mut app = App {
            sim: Simulation::new(nbody_core::system::NBodySystem::si(), Default::default()),
            scenario: scenario.clone(),
            scenario_name: String::new(),
            units: UnitSystem::SI,
            unit_spec: UnitSpec::default(),
            projectile_template: None,
            running: false,
            time_scale: 1.0,
            max_steps_per_frame: if cfg!(target_arch = "wasm32") { 6_000 } else { 20_000 },
            steps_last_frame: 0,
            lagging: false,
            view: View { center: Vec3::ZERO, scale: 1.0, yaw: 0.0, pitch: 0.0 },
            follow: None,
            show_labels: true,
            show_trails: true,
            trail_len: 2500,
            trails: HashMap::new(),
            selected: None,
            edit_loaded_for: None,
            edit_name: String::new(),
            shell_thickness: 0.0,
            shell_material: Material::Ice,
            orbit_target: None,
            new_name: "New planet".into(),
            new_mass: 1.0,
            new_radius: 0.0,
            new_material: Material::Rock,
            new_core: Some(Material::Iron),
            new_kind: BodyKind::Planet,
            new_around: None,
            new_a: 1.0,
            new_e: 0.0,
            new_i: 0.0,
            new_pos: [0.0; 3],
            new_vel: [0.0; 3],
            p_name: "Projectile".into(),
            p_v_inf: 1.0,
            p_b: 1.0,
            p_mass: 1.0,
            p_material: Material::Ice,
            p_angle: 0.0,
            p_side: true,
            p_start: 0.0,
            p_max_time: 100.0,
            sweep_b: [0.5, 5.0],
            sweep_n: 24,
            sweep_v: [0.5, 2.0],
            sweep_nv: 4,
            sweep_result: String::new(),
            fixed_dt: 1.0,
            adaptive: false,
            eta: 0.01,
            preset_index: 0,
            show_bodies: true,
            show_projectile: true,
            layout_initialised: false,
            file_path: "scenario.toml".into(),
            energy_history: VecDeque::new(),
            log: VecDeque::new(),
            error: None,
        };
        app.load_scenario(scenario, "Solar System".into());
        app
    }

    // ---- scenario handling ------------------------------------------------------------

    fn load_scenario(&mut self, scenario: Scenario, name: String) {
        match scenario.build() {
            Ok(built) => {
                self.units = built.units;
                self.unit_spec = built.unit_spec.clone();
                self.projectile_template = built.projectile.clone();
                match built.config.timestep {
                    TimestepPolicy::Fixed { dt } => {
                        self.fixed_dt = dt / self.units.time;
                        self.adaptive = false;
                    }
                    TimestepPolicy::Adaptive { eta, dt_max, .. } => {
                        self.fixed_dt = if dt_max.is_finite() { dt_max / self.units.time } else { 1.0 };
                        self.adaptive = true;
                        self.eta = eta;
                    }
                }
                let duration = built.duration.unwrap_or(10.0 * self.units.time);
                self.time_scale = (duration / self.units.time / 30.0).max(1e-6);
                if let Some(p) = &built.projectile {
                    self.p_name = p.name.clone();
                    self.p_v_inf = p.v_infinity / self.units.velocity();
                    self.p_b = p.impact_parameter / self.units.length;
                    self.p_mass = p.structure.mass() / self.units.mass;
                    self.p_material = p.structure.dominant_material().clone();
                    self.p_angle = p.direction.y.atan2(p.direction.x).to_degrees();
                    self.p_side = p.offset_direction.map(|o| Vec3::Z.cross(p.direction).dot(o) >= 0.0).unwrap_or(true);
                    self.p_start = p.start_distance.map(|d| d / self.units.length).unwrap_or(0.0);
                    self.p_max_time = built.projectile_max_time.map(|t| t / self.units.time).unwrap_or(100.0);
                    self.sweep_b = [0.5 * self.p_b, 1.5 * self.p_b];
                    self.sweep_v = [0.5 * self.p_v_inf, 2.0 * self.p_v_inf];
                } else {
                    let ext = built.system.extent() / self.units.length;
                    self.p_b = (0.5 * ext).max(1e-3);
                    self.p_v_inf = 1.0;
                    self.p_mass = built.system.bodies().iter().map(|b| b.mass()).fold(f64::INFINITY, f64::min) / self.units.mass;
                    self.p_start = 0.0;
                    self.sweep_b = [0.1 * ext, 2.0 * ext];
                    self.sweep_v = [0.5, 2.0];
                }
                self.new_mass = built.system.bodies().iter().map(|b| b.mass()).fold(f64::INFINITY, f64::min) / self.units.mass;
                self.new_a = built.system.extent() / self.units.length * 0.5;
                self.new_around = built.system.most_massive().map(|b| b.id());
                self.orbit_target = self.new_around;
                self.sim = Simulation::new(built.system, built.config);
                self.scenario = scenario;
                self.scenario_name = name;
                self.running = false;
                self.selected = None;
                self.edit_loaded_for = None;
                self.follow = None;
                self.trails.clear();
                self.energy_history.clear();
                self.sweep_result.clear();
                self.error = None;
                self.log_line(format!("loaded '{}' ({} bodies)", self.scenario_name, self.sim.system.len()));
                self.fit_view();
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    fn reset(&mut self) {
        let sc = self.scenario.clone();
        let name = self.scenario_name.clone();
        self.load_scenario(sc, name);
    }

    fn log_line(&mut self, s: String) {
        self.log.push_back(s);
        while self.log.len() > 200 {
            self.log.pop_front();
        }
    }

    fn fit_view(&mut self) {
        let ext = self.sim.system.extent().max(self.sim.system.bodies().iter().map(|b| b.radius()).fold(0.0, f64::max) * 4.0);
        self.view.center = self.sim.system.center_of_mass();
        self.view.scale = if ext > 0.0 { 380.0 / ext } else { 1.0 };
    }

    // ---- unit helpers -----------------------------------------------------------------

    fn fmt_len(&self, m: f64) -> String {
        format!("{} {}", sig(m / self.units.length), self.unit_spec.length)
    }
    fn fmt_mass(&self, kg: f64) -> String {
        format!("{} {}", sig(kg / self.units.mass), self.unit_spec.mass)
    }
    fn fmt_time(&self, s: f64) -> String {
        format!("{} {}", sig(s / self.units.time), self.unit_spec.time)
    }
    fn fmt_speed(&self, v: f64) -> String {
        format!("{} {}/{}", sig(v / self.units.velocity()), self.unit_spec.length, self.unit_spec.time)
    }

    // ---- stepping ---------------------------------------------------------------------

    fn advance_frame(&mut self, real_dt: f64) {
        let budget = self.time_scale * self.units.time * real_dt.clamp(0.0, 0.1);
        let target = self.sim.system.time + budget;
        let nominal = self.sim.config().timestep.nominal_dt().max(1e-300);
        let expected = (budget / nominal).ceil() as usize;
        let stride = (expected / 60).max(1);
        let mut steps = 0;
        while self.sim.system.time < target && steps < self.max_steps_per_frame {
            let dt = self.sim.choose_dt().min(target - self.sim.system.time);
            if dt <= 0.0 {
                break;
            }
            let n_before = self.sim.collisions.len();
            self.sim.step_with(dt);
            if self.sim.collisions.len() > n_before {
                for e in self.sim.collisions[n_before..].to_vec() {
                    self.log_line(format!("t = {}: {}", self.fmt_time(e.time), e.describe()));
                }
            }
            steps += 1;
            if steps % stride == 0 {
                self.push_trails();
            }
        }
        self.lagging = steps >= self.max_steps_per_frame;
        self.steps_last_frame = steps;
        let d = self.sim.diagnostics();
        let e = self.sim.conservation.record(&d) as f32;
        self.energy_history.push_back(e);
        while self.energy_history.len() > 600 {
            self.energy_history.pop_front();
        }
    }

    fn push_trails(&mut self) {
        if !self.show_trails {
            return;
        }
        let ids: Vec<BodyId> = self.sim.system.ids();
        self.trails.retain(|id, _| ids.contains(id));
        for b in self.sim.system.bodies() {
            let t = self.trails.entry(b.id()).or_default();
            t.push_back(b.position);
            while t.len() > self.trail_len {
                t.pop_front();
            }
        }
    }

    fn apply_timestep(&mut self) {
        let policy = if self.adaptive {
            TimestepPolicy::Adaptive { eta: self.eta, dt_min: 0.0, dt_max: self.fixed_dt * self.units.time }
        } else {
            TimestepPolicy::Fixed { dt: self.fixed_dt * self.units.time }
        };
        self.sim.set_timestep(policy);
    }

    // ---- panels -----------------------------------------------------------------------

    fn top_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            let names = presets::names();
            egui::ComboBox::from_id_salt("preset")
                .selected_text(names[self.preset_index].0)
                .width(170.0)
                .show_ui(ui, |ui| {
                    for (i, (n, d)) in names.iter().enumerate() {
                        ui.selectable_value(&mut self.preset_index, i, *n).on_hover_text(*d);
                    }
                });
            if ui.button("Load preset").clicked() {
                let (n, _) = names[self.preset_index];
                if let Some(sc) = presets::preset(n) {
                    let name = sc.name.clone();
                    self.load_scenario(sc, name);
                }
            }
            #[cfg(not(target_arch = "wasm32"))]
            {
            ui.add(egui::TextEdit::singleline(&mut self.file_path).desired_width(160.0));
            if ui.button("Open file").clicked() {
                match Scenario::load(&self.file_path) {
                    Ok(sc) => {
                        let name = if sc.name.is_empty() { self.file_path.clone() } else { sc.name.clone() };
                        self.load_scenario(sc, name);
                    }
                    Err(e) => self.error = Some(e.to_string()),
                }
            }
            if ui.button("Save file").clicked() {
                match Scenario::snapshot(&self.sim.system, self.sim.config(), &self.unit_spec, &self.scenario_name).and_then(|s| s.save(&self.file_path)) {
                    Ok(()) => self.log_line(format!("saved {}", self.file_path)),
                    Err(e) => self.error = Some(e.to_string()),
                }
            }
            }
            ui.separator();
            if ui.button(if self.running { "⏸ Pause (space)" } else { "▶ Run (space)" }).clicked() {
                self.running = !self.running;
            }
            if ui.button("Step").clicked() {
                self.sim.step();
                self.push_trails();
            }
            if ui.button("Reset").clicked() {
                self.reset();
            }
            ui.separator();
            ui.label("speed");
            ui.add(
                egui::Slider::new(&mut self.time_scale, 1e-4..=1e5)
                    .logarithmic(true)
                    .suffix(format!(" {}/s", self.unit_spec.time)),
            );
            ui.separator();
            let mut integ = self.sim.config().integrator;
            egui::ComboBox::from_id_salt("integrator").selected_text(integ.name()).show_ui(ui, |ui| {
                for k in IntegratorKind::ALL {
                    ui.selectable_value(&mut integ, k, k.name()).on_hover_text(k.description());
                }
            });
            if integ != self.sim.config().integrator {
                self.sim.set_integrator(integ);
            }
            let mut coll = self.sim.config().collisions;
            egui::ComboBox::from_id_salt("collisions").selected_text(format!("collisions: {coll}")).show_ui(ui, |ui| {
                for m in CollisionMode::ALL {
                    ui.selectable_value(&mut coll, m, m.name());
                }
            });
            if coll != self.sim.config().collisions {
                self.sim.set_collisions(coll);
            }
            ui.separator();
            ui.label(format!("dt [{}]", self.unit_spec.time));
            let changed = drag(ui, &mut self.fixed_dt, 1e-12).changed();
            let toggled = ui.checkbox(&mut self.adaptive, "adaptive").on_hover_text("dt = η · (shortest two-body dynamical time), capped at the value on the left").changed();
            let eta_changed = self.adaptive && ui.add(egui::DragValue::new(&mut self.eta).speed(0.001).range(1e-4..=1.0).prefix("η = ")).changed();
            if (changed || toggled || eta_changed) && self.fixed_dt > 0.0 {
                self.apply_timestep();
            }
            ui.separator();
            if ui.button("Barycentre").clicked() {
                self.sim.system.move_to_barycentre();
                self.sim.reset_baseline();
            }
            if ui.button("Fit view (F)").clicked() {
                self.fit_view();
            }
            ui.checkbox(&mut self.show_labels, "labels");
            ui.checkbox(&mut self.show_trails, "trails");
            ui.separator();
            ui.checkbox(&mut self.show_bodies, "bodies panel");
            ui.checkbox(&mut self.show_projectile, "projectile panel");
        });
    }

    fn bodies_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading(format!("Bodies ({})", self.sim.system.len()));
        let rows: Vec<(BodyId, String, [u8; 3], f64, &'static str)> =
            self.sim.system.bodies().iter().map(|b| (b.id(), b.name.clone(), b.colour(), b.mass(), b.kind.name())).collect();
        egui::ScrollArea::vertical().id_salt("bodylist").max_height(200.0).show(ui, |ui| {
            for (id, name, col, mass, kind) in &rows {
                ui.horizontal(|ui| {
                    let (r, _) = ui.allocate_exact_size(Vec2::splat(12.0), Sense::hover());
                    ui.painter().circle_filled(r.center(), 5.0, c32(*col));
                    let label = format!("{name}  ({kind}, {})", self.fmt_mass(*mass));
                    if ui.selectable_label(self.selected == Some(*id), label).clicked() {
                        self.selected = Some(*id);
                    }
                });
            }
        });
        ui.separator();
        if let Some(id) = self.selected {
            if self.sim.system.get(id).is_none() {
                self.selected = None;
            }
        }
        if let Some(id) = self.selected {
            self.body_editor(ui, id);
        } else {
            ui.label("Click a body in the list or on the canvas to edit it.");
        }
        ui.separator();
        ui.collapsing("Add a body", |ui| self.add_body_ui(ui));
    }

    fn material_combo(ui: &mut egui::Ui, id: impl std::hash::Hash, current: &mut Material) -> bool {
        let before = current.clone();
        egui::ComboBox::from_id_salt(id).selected_text(current.name().to_string()).show_ui(ui, |ui| {
            for m in Material::BUILTIN.iter() {
                ui.selectable_value(current, m.clone(), m.name()).on_hover_text(format!("{} ({} kg/m³, e = {:.2})", m.description(), m.density(), m.restitution()));
            }
        });
        *current != before
    }

    fn kind_combo(ui: &mut egui::Ui, id: impl std::hash::Hash, current: &mut BodyKind) -> bool {
        let before = *current;
        egui::ComboBox::from_id_salt(id).selected_text(current.name()).show_ui(ui, |ui| {
            for k in [BodyKind::Star, BodyKind::Planet, BodyKind::Moon, BodyKind::Asteroid, BodyKind::Comet, BodyKind::Projectile, BodyKind::Remnant, BodyKind::Particle] {
                ui.selectable_value(current, k, k.name());
            }
        });
        *current != before
    }

    fn body_combo(&self, ui: &mut egui::Ui, id: &str, current: &mut Option<BodyId>, exclude: Option<BodyId>) {
        let name_of = |b: Option<BodyId>| b.and_then(|i| self.sim.system.get(i)).map(|b| b.name.clone()).unwrap_or_else(|| "(none)".into());
        egui::ComboBox::from_id_salt(id).selected_text(name_of(*current)).show_ui(ui, |ui| {
            for b in self.sim.system.bodies() {
                if Some(b.id()) == exclude {
                    continue;
                }
                ui.selectable_value(current, Some(b.id()), b.name.clone());
            }
        });
    }

    fn body_editor(&mut self, ui: &mut egui::Ui, id: BodyId) {
        let body = self.sim.system.get(id).unwrap().clone();
        if self.edit_loaded_for != Some(id) {
            self.edit_name = body.name.clone();
            self.edit_loaded_for = Some(id);
        }
        ui.heading(&body.name);
        let (u_len, u_mass, u_vel) = (self.units.length, self.units.mass, self.units.velocity());
        let mut changed = false;
        let mut new_mass = body.mass() / u_mass;
        let mut new_radius = body.radius() / u_len;
        let mut new_material = body.material().clone();
        let mut new_kind = body.kind;
        let mut pos = [body.position.x / u_len, body.position.y / u_len, body.position.z / u_len];
        let mut vel = [body.velocity.x / u_vel, body.velocity.y / u_vel, body.velocity.z / u_vel];
        let (mut mass_changed, mut radius_changed, mut material_changed, mut pos_changed, mut vel_changed, mut rename) = (false, false, false, false, false, false);
        egui::Grid::new(("edit", id)).num_columns(2).spacing([8.0, 4.0]).show(ui, |ui| {
            ui.label("name");
            if ui.text_edit_singleline(&mut self.edit_name).lost_focus() {
                rename = true;
            }
            ui.end_row();
            ui.label(format!("mass [{}]", self.unit_spec.mass));
            mass_changed = drag(ui, &mut new_mass, 1e-30).on_hover_text("homologous rescaling: densities fixed, radius ∝ M^(1/3)").changed();
            ui.end_row();
            ui.label(format!("radius [{}]", self.unit_spec.length));
            radius_changed = drag(ui, &mut new_radius, 1e-30).on_hover_text("homologous rescaling: densities fixed, mass ∝ R³").changed();
            ui.end_row();
            ui.label("material");
            material_changed = Self::material_combo(ui, ("mat", id), &mut new_material);
            ui.end_row();
            ui.label("kind");
            if Self::kind_combo(ui, ("kind", id), &mut new_kind) {
                changed = true;
            }
            ui.end_row();
            ui.label(format!("position [{}]", self.unit_spec.length));
            ui.horizontal(|ui| {
                for p in &mut pos {
                    pos_changed |= drag(ui, p, 1e-6).changed();
                }
            });
            ui.end_row();
            ui.label(format!("velocity [{}/{}]", self.unit_spec.length, self.unit_spec.time));
            ui.horizontal(|ui| {
                for v in &mut vel {
                    vel_changed |= drag(ui, v, 1e-6).changed();
                }
            });
            ui.end_row();
        });
        {
            let g = self.sim.system.g;
            let b = self.sim.system.get_mut(id).unwrap();
            if rename && !self.edit_name.trim().is_empty() {
                b.name = self.edit_name.trim().to_string();
                changed = true;
            }
            if mass_changed && new_mass > 0.0 {
                let _ = b.set_mass(new_mass * u_mass);
                changed = true;
            }
            if radius_changed && new_radius > 0.0 {
                let _ = b.set_radius(new_radius * u_len);
                changed = true;
            }
            if material_changed {
                let _ = b.set_material(new_material);
                changed = true;
            }
            if new_kind != b.kind {
                b.kind = new_kind;
            }
            if pos_changed {
                b.position = Vec3::new(pos[0], pos[1], pos[2]) * u_len;
                changed = true;
            }
            if vel_changed {
                b.velocity = Vec3::new(vel[0], vel[1], vel[2]) * u_vel;
                changed = true;
            }
            let _ = g;
        }
        let body = self.sim.system.get(id).unwrap().clone();
        let s = body.structure();
        ui.label(format!(
            "ρ̄ = {:.0} kg/m³   I/MR² = {:.3}   v_esc = {}   g = {:.3} m/s²",
            s.mean_density(),
            s.inertia_factor(),
            self.fmt_speed(body.escape_velocity(self.sim.system.g)),
            body.surface_gravity(self.sim.system.g)
        ));
        ui.collapsing(format!("Internal structure ({} layer{})", s.layers().len(), if s.layers().len() == 1 { "" } else { "s" }), |ui| {
            let total = s.mass();
            let mut inner = 0.0;
            for l in s.layers() {
                let m = l.material.density() * (structure::sphere_volume(l.outer_radius) - structure::sphere_volume(inner));
                ui.horizontal(|ui| {
                    let (r, _) = ui.allocate_exact_size(Vec2::splat(12.0), Sense::hover());
                    ui.painter().circle_filled(r.center(), 5.0, c32(l.material.colour()));
                    ui.label(format!("r ≤ {}  {}  {:.0} kg/m³  {:.1}% of mass", self.fmt_len(l.outer_radius), l.material.name(), l.material.density(), 100.0 * m / total));
                });
                inner = l.outer_radius;
            }
            ui.horizontal(|ui| {
                ui.label("add shell: thickness");
                if self.shell_thickness <= 0.0 {
                    self.shell_thickness = 0.1 * body.radius() / u_len;
                }
                drag(ui, &mut self.shell_thickness, 1e-30);
                Self::material_combo(ui, ("shellmat", id), &mut self.shell_material);
                if ui.button("Add").clicked() && self.shell_thickness > 0.0 {
                    let mut ns = s.clone();
                    if ns.push_shell(self.shell_thickness * u_len, self.shell_material.clone()).is_ok() {
                        self.sim.system.get_mut(id).unwrap().set_structure(ns);
                        changed = true;
                    }
                }
                if ui.button("Remove outer shell").clicked() {
                    let mut ns = s.clone();
                    if ns.pop_shell().is_some() {
                        self.sim.system.get_mut(id).unwrap().set_structure(ns);
                        changed = true;
                    }
                }
            });
        });
        ui.horizontal(|ui| {
            ui.label("orbit around");
            let mut target = self.orbit_target;
            self.body_combo(ui, "orbit_target", &mut target, Some(id));
            self.orbit_target = target;
            if let Some(t) = target {
                if ui.button("Circularise").on_hover_text("set the velocity for a circular orbit at the current distance").clicked() {
                    if self.sim.system.set_circular_orbit(id, t, Vec3::Z).is_ok() {
                        changed = true;
                    }
                }
            }
            let following = self.follow == Some(id);
            if ui.selectable_label(following, "Follow").clicked() {
                self.follow = if following { None } else { Some(id) };
            }
            if ui.button("Delete").clicked() {
                self.sim.system.remove(id);
                self.selected = None;
                changed = true;
            }
        });
        if let Some(t) = self.orbit_target {
            if let (Some(p), Some(b)) = (self.sim.system.get(t), self.sim.system.get(id)) {
                if p.id() != b.id() {
                    let mu = self.sim.system.g * (p.mass() + b.mass());
                    let el = orbit::elements_from_state(mu, b.position - p.position, b.velocity - p.velocity);
                    let hill = orbit::hill_radius(el.elements.a.abs(), el.elements.e.min(0.999), b.mass(), p.mass());
                    let roche = orbit::roche_limit_rigid(p.radius(), p.mean_density(), b.mean_density());
                    let txt = if el.is_bound() {
                        format!("a = {}, e = {:.4}, i = {:.2}°, T = {}", self.fmt_len(el.elements.a), el.elements.e, el.elements.inclination.to_degrees(), self.fmt_time(el.period.unwrap_or(f64::NAN)))
                    } else {
                        format!("unbound: e = {:.3}, v∞ = {}", el.elements.e, self.fmt_speed(el.v_infinity.unwrap_or(0.0)))
                    };
                    ui.small(format!("relative to {}: {}", p.name, txt));
                    ui.small(format!("Hill radius {}, rigid Roche limit {}", self.fmt_len(hill), self.fmt_len(roche)));
                }
            }
        }
        if changed {
            self.sim.reset_baseline();
            self.trails.remove(&id);
        }
    }

    fn add_body_ui(&mut self, ui: &mut egui::Ui) {
        egui::Grid::new("addbody").num_columns(2).spacing([8.0, 4.0]).show(ui, |ui| {
            ui.label("name");
            ui.text_edit_singleline(&mut self.new_name);
            ui.end_row();
            ui.label(format!("mass [{}]", self.unit_spec.mass));
            drag(ui, &mut self.new_mass, 1e-30);
            ui.end_row();
            ui.label(format!("radius [{}] (0 = from density)", self.unit_spec.length));
            drag(ui, &mut self.new_radius, 1e-9);
            ui.end_row();
            ui.label("material");
            Self::material_combo(ui, "newmat", &mut self.new_material);
            ui.end_row();
            ui.label("core (with radius)");
            ui.horizontal(|ui| {
                let mut has_core = self.new_core.is_some();
                if ui.checkbox(&mut has_core, "").changed() {
                    self.new_core = if has_core { Some(Material::Iron) } else { None };
                }
                if let Some(c) = &mut self.new_core {
                    Self::material_combo(ui, "newcore", c);
                }
            });
            ui.end_row();
            ui.label("kind");
            Self::kind_combo(ui, "newkind", &mut self.new_kind);
            ui.end_row();
            ui.label("orbit around");
            let mut around = self.new_around;
            self.body_combo(ui, "newaround", &mut around, None);
            self.new_around = around;
            ui.end_row();
            if self.new_around.is_some() {
                ui.label(format!("a [{}]", self.unit_spec.length));
                drag(ui, &mut self.new_a, 1e-9);
                ui.end_row();
                ui.label("e");
                ui.add(egui::DragValue::new(&mut self.new_e).speed(0.005).range(0.0..=0.999));
                ui.end_row();
                ui.label("inclination [°]");
                ui.add(egui::DragValue::new(&mut self.new_i).speed(0.5));
                ui.end_row();
            } else {
                ui.label(format!("position [{}]", self.unit_spec.length));
                ui.horizontal(|ui| {
                    for p in &mut self.new_pos {
                        drag(ui, p, 1e-6);
                    }
                });
                ui.end_row();
                ui.label("velocity");
                ui.horizontal(|ui| {
                    for v in &mut self.new_vel {
                        drag(ui, v, 1e-6);
                    }
                });
                ui.end_row();
            }
        });
        ui.horizontal(|ui| {
            if ui.button("Add body").clicked() {
                match self.build_new_body() {
                    Ok(id) => {
                        self.selected = Some(id);
                        self.sim.reset_baseline();
                        self.log_line(format!("added {}", self.new_name));
                    }
                    Err(e) => self.error = Some(e),
                }
            }
            if ui.button("Clear all").clicked() {
                self.sim.system.clear();
                self.trails.clear();
                self.selected = None;
                self.sim.reset_baseline();
            }
        });
    }

    fn build_new_body(&mut self) -> Result<BodyId, String> {
        let u = self.units;
        let mass = self.new_mass * u.mass;
        let structure = if self.new_radius > 0.0 {
            let r = self.new_radius * u.length;
            match &self.new_core {
                Some(core) => structure::fit_two_layer(core.clone(), self.new_material.clone(), mass, r),
                None => {
                    let density = mass / structure::sphere_volume(r);
                    SphericalStructure::uniform_with_radius(Material::custom(self.new_material.name(), density, self.new_material.restitution()), r)
                }
            }
        } else {
            SphericalStructure::uniform_with_mass(self.new_material.clone(), mass)
        }
        .map_err(|e| e.to_string())?;
        let name = if self.new_name.trim().is_empty() { format!("body{}", self.sim.system.len() + 1) } else { self.new_name.trim().to_string() };
        let mut body = Body::new(name, structure, Vec3::from(self.new_pos) * u.length, Vec3::from(self.new_vel) * u.velocity()).with_kind(self.new_kind);
        if let Some(around) = self.new_around {
            let primary = self.sim.system.get(around).ok_or("orbit target no longer exists")?.clone();
            let mu = self.sim.system.g * (primary.mass() + body.mass());
            let el = Elements { a: self.new_a * u.length, e: self.new_e, inclination: self.new_i.to_radians(), raan: 0.0, arg_periapsis: 0.0, true_anomaly: 0.0 };
            let (r, v) = orbit::state_from_elements(mu, &el).map_err(|e| e.to_string())?;
            body.position = primary.position + r;
            body.velocity = primary.velocity + v;
        }
        Ok(self.sim.system.add(body))
    }

    fn projectile_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Projectile mode");
        ui.small("An interloper arrives from infinity on the exact two-body hyperbola set by v∞ and the impact parameter b. The tracker then decides whether the system captures it, ejects it, or swallows it, and whether it kicked anyone else out.");
        let (ul, ut, um) = (self.unit_spec.length.clone(), self.unit_spec.time.clone(), self.unit_spec.mass.clone());
        egui::Grid::new("proj").num_columns(2).spacing([8.0, 4.0]).show(ui, |ui| {
            ui.label("name");
            ui.text_edit_singleline(&mut self.p_name);
            ui.end_row();
            ui.label(format!("v∞ [{ul}/{ut}]"));
            drag(ui, &mut self.p_v_inf, 1e-6);
            ui.end_row();
            ui.label(format!("impact parameter b [{ul}]"));
            drag(ui, &mut self.p_b, 1e-9);
            ui.end_row();
            ui.label(format!("mass [{um}]"));
            drag(ui, &mut self.p_mass, 1e-30);
            ui.end_row();
            ui.label("material");
            Self::material_combo(ui, "pmat", &mut self.p_material);
            ui.end_row();
            ui.label("approach direction [°]");
            ui.add(egui::DragValue::new(&mut self.p_angle).speed(1.0)).on_hover_text("direction of the velocity at infinity, in the x-y plane");
            ui.end_row();
            ui.label("pass on the");
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.p_side, true, "left (prograde)");
                ui.selectable_value(&mut self.p_side, false, "right (retrograde)");
            });
            ui.end_row();
            ui.label(format!("start distance [{ul}] (0 = auto)"));
            drag(ui, &mut self.p_start, 1e-6);
            ui.end_row();
            ui.label(format!("give up after [{ut}]"));
            drag(ui, &mut self.p_max_time, 1e-6);
            ui.end_row();
        });
        if !self.sim.system.is_empty() {
            let mu = self.sim.system.g * (self.sim.system.total_mass() + self.p_mass * self.units.mass);
            let v = self.p_v_inf * self.units.velocity();
            let b = self.p_b * self.units.length;
            ui.small(format!(
                "two-body prediction about the barycentre: periapsis {}, deflection {:.1}°",
                self.fmt_len(orbit::hyperbolic_periapsis(mu, v, b)),
                orbit::hyperbolic_deflection(mu, v, b.max(1e-300)).to_degrees()
            ));
        }
        ui.horizontal(|ui| {
            if ui.button("🚀 Launch").clicked() {
                match self.projectile_spec() {
                    Ok(spec) => match self.sim.inject_projectile(&spec) {
                        Ok(id) => {
                            self.selected = Some(id);
                            self.running = true;
                            self.log_line(format!("launched {} (v∞ = {}, b = {})", spec.name, self.fmt_speed(spec.v_infinity), self.fmt_len(spec.impact_parameter)));
                            self.fit_view();
                        }
                        Err(e) => self.error = Some(e.to_string()),
                    },
                    Err(e) => self.error = Some(e),
                }
            }
            if ui.button("Remove projectiles").clicked() {
                let ids: Vec<BodyId> = self.sim.trackers.iter().map(|t| t.id).collect();
                for id in ids {
                    self.sim.system.remove(id);
                }
                self.sim.trackers.clear();
                self.sim.reset_baseline();
            }
        });
        if let Some(t) = self.sim.trackers.last().cloned() {
            ui.separator();
            let (label, colour) = match &t.outcome {
                Outcome::Captured { .. } => ("CAPTURED", Color32::from_rgb(90, 220, 120)),
                Outcome::Ejected { .. } => ("EJECTED", Color32::from_rgb(240, 170, 60)),
                Outcome::Merged { .. } => ("MERGED", Color32::from_rgb(230, 80, 80)),
                Outcome::InFlight { phase } => (
                    match phase {
                        Phase::Approaching => "approaching",
                        Phase::Receding => "receding, unbound",
                        Phase::Interacting => "falling back in",
                        Phase::Encounter => "close encounter",
                    },
                    Color32::LIGHT_GRAY,
                ),
                Outcome::Undetermined => ("undetermined", Color32::GRAY),
            };
            ui.label(egui::RichText::new(format!("{}: {label}", t.name)).strong().color(colour).size(16.0));
            match &t.outcome {
                Outcome::Captured { orbits, elements } => {
                    ui.label(format!("a = {}, e = {:.4}, period {}, {} full orbits so far", self.fmt_len(elements.elements.a), elements.elements.e, self.fmt_time(elements.period.unwrap_or(f64::NAN)), orbits));
                }
                Outcome::Ejected { v_infinity_out, deflection, specific_energy_change } => {
                    ui.label(format!("leaves at v∞ = {} (arrived at {}), deflected {:.1}°, Δε = {:+.3e} J/kg", self.fmt_speed(*v_infinity_out), self.fmt_speed(t.incoming_v_infinity), deflection.to_degrees(), specific_energy_change));
                }
                Outcome::Merged { into, time } => {
                    ui.label(format!("absorbed by {into} at t = {}", self.fmt_time(*time)));
                }
                _ => {}
            }
            ui.small(format!("closest approach to barycentre {} at t = {}", self.fmt_len(t.closest_approach), self.fmt_time(t.closest_approach_time)));
            if let Some((name, gap)) = &t.closest_body {
                ui.small(format!("closest surface gap {} to {name}", self.fmt_len(*gap)));
            }
            ui.small(format!("specific energy {:+.3e} J/kg (was {:+.3e}); periapsis passages {}, apoapsis passages {}", t.specific_energy, t.initial_specific_energy, t.periapsis_passages, t.apoapsis_passages));
            if !t.displaced.is_empty() {
                ui.label(egui::RichText::new(format!("unbound by the encounter: {}", t.displaced.join(", "))).color(Color32::from_rgb(240, 120, 120)));
            }
            if !t.absorbed.is_empty() {
                ui.small(format!("absorbed: {}", t.absorbed.join(", ")));
            }
        }
        ui.separator();
        ui.collapsing("Outcome sweep", |ui| {
            ui.small("Runs one experiment per (v∞, b) on a copy of the current system. C = captured, . = ejected, # = merged, ? = undetermined.");
            egui::Grid::new("sweep").num_columns(2).show(ui, |ui| {
                ui.label(format!("b range [{ul}]"));
                ui.horizontal(|ui| {
                    drag(ui, &mut self.sweep_b[0], 1e-9);
                    drag(ui, &mut self.sweep_b[1], 1e-9);
                    ui.add(egui::DragValue::new(&mut self.sweep_n).range(1..=200));
                });
                ui.end_row();
                ui.label(format!("v∞ range [{ul}/{ut}]"));
                ui.horizontal(|ui| {
                    drag(ui, &mut self.sweep_v[0], 1e-9);
                    drag(ui, &mut self.sweep_v[1], 1e-9);
                    ui.add(egui::DragValue::new(&mut self.sweep_nv).range(1..=20));
                });
                ui.end_row();
            });
            if ui.button("Run sweep (blocks the UI)").clicked() {
                self.run_sweep();
            }
            if !self.sweep_result.is_empty() {
                ui.monospace(&self.sweep_result);
            }
        });
    }

    fn projectile_spec(&self) -> Result<ProjectileSpec, String> {
        let u = self.units;
        if self.p_mass <= 0.0 {
            return Err("projectile mass must be positive".into());
        }
        let structure = SphericalStructure::uniform_with_mass(self.p_material.clone(), self.p_mass * u.mass).map_err(|e| e.to_string())?;
        let (s, c) = self.p_angle.to_radians().sin_cos();
        let dir = Vec3::new(c, s, 0.0);
        let left = Vec3::Z.cross(dir);
        let mut spec = ProjectileSpec::new(self.p_name.clone(), structure, self.p_v_inf * u.velocity(), self.p_b * u.length)
            .with_direction(dir)
            .with_offset_direction(if self.p_side { left } else { -left });
        if self.p_start > 0.0 {
            spec = spec.with_start_distance(self.p_start * u.length);
        }
        Ok(spec)
    }

    fn run_sweep(&mut self) {
        let template = match self.projectile_spec() {
            Ok(s) => s,
            Err(e) => {
                self.error = Some(e);
                return;
            }
        };
        let u = self.units;
        let n = self.sweep_n.max(1);
        let nv = self.sweep_nv.max(1);
        let bs: Vec<f64> = (0..n).map(|k| (self.sweep_b[0] + (self.sweep_b[1] - self.sweep_b[0]) * k as f64 / (n.max(2) - 1) as f64) * u.length).collect();
        let vs: Vec<f64> = (0..nv).map(|k| (self.sweep_v[0] + (self.sweep_v[1] - self.sweep_v[0]) * k as f64 / (nv.max(2) - 1) as f64) * u.velocity()).collect();
        let r = projectile::sweep(&self.sim.system, self.sim.config(), &template, &bs, &vs, self.p_max_time * u.time, None);
        let fr = r.cross_section_fractions();
        let mut out = format!("b from {} to {} ({} values)\n", sig(bs[0] / u.length), sig(bs[n - 1] / u.length), n);
        for (iv, row) in r.grid.iter().enumerate() {
            let line: String = row.iter().map(|c| c.symbol()).collect();
            out.push_str(&format!("v∞={:>8} | {line} | capture {:>4.0}%  merge {:>4.0}%\n", sig(vs[iv] / u.velocity()), 100.0 * fr[iv][0], 100.0 * fr[iv][2]));
        }
        self.sweep_result = out;
    }

    fn bottom_panel(&mut self, ui: &mut egui::Ui) {
        let d = self.sim.diagnostics();
        let e = self.sim.conservation.energy_error(&d);
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(format!("t = {}", self.fmt_time(d.time))).strong());
            ui.separator();
            ui.label(format!("{} bodies, {} steps, {} force evals", d.bodies, self.sim.steps, self.sim.force_evaluations));
            ui.separator();
            ui.label(format!("last dt = {}", self.fmt_time(self.sim.last_dt)));
            ui.separator();
            ui.label(format!("E = {:+.4e} J, 2K/|U| = {:.3}", d.total, d.virial_ratio));
            ui.separator();
            let colour = if e.abs() < 1e-6 { Color32::from_rgb(90, 220, 120) } else if e.abs() < 1e-3 { Color32::from_rgb(240, 200, 60) } else { Color32::from_rgb(230, 80, 80) };
            ui.label(egui::RichText::new(format!("ΔE/E = {e:+.2e} (max {:.1e})", self.sim.conservation.max_energy_error)).color(colour));
            ui.label(format!("Δp = {:.1e}, ΔL = {:.1e}", self.sim.conservation.momentum_error(&d), self.sim.conservation.angular_momentum_error(&d)));
            if self.lagging {
                ui.label(egui::RichText::new(format!("⚠ capped at {} steps/frame: lower the speed or raise dt", self.max_steps_per_frame)).color(Color32::from_rgb(240, 170, 60)));
            }
        });
        ui.horizontal(|ui| {
            let (rect, painter) = ui.allocate_painter(Vec2::new(300.0, 46.0), Sense::hover());
            let rect = rect.rect;
            painter.rect_filled(rect, 3.0f32, Color32::from_rgb(25, 25, 30));
            if self.energy_history.len() > 1 {
                let max = self.energy_history.iter().fold(1e-16f32, |m, v| m.max(v.abs()));
                let n = self.energy_history.len();
                let pts: Vec<Pos2> = self
                    .energy_history
                    .iter()
                    .enumerate()
                    .map(|(i, v)| Pos2::new(rect.left() + rect.width() * i as f32 / (n - 1) as f32, rect.center().y - 0.45 * rect.height() * v / max))
                    .collect();
                painter.line_segment([Pos2::new(rect.left(), rect.center().y), Pos2::new(rect.right(), rect.center().y)], Stroke::new(1.0f32, Color32::from_gray(60)));
                painter.add(egui::Shape::line(pts, Stroke::new(1.0f32, Color32::from_rgb(120, 200, 255))));
                painter.text(rect.left_top() + Vec2::new(4.0, 2.0), Align2::LEFT_TOP, format!("ΔE/E, ±{max:.1e}"), FontId::proportional(10.0), Color32::GRAY);
            }
            ui.vertical(|ui| {
                egui::ScrollArea::vertical().id_salt("log").max_height(46.0).stick_to_bottom(true).show(ui, |ui| {
                    for line in self.log.iter().rev().take(6).collect::<Vec<_>>().into_iter().rev() {
                        ui.small(line);
                    }
                });
            });
        });
        if let Some(err) = self.error.clone() {
            ui.horizontal(|ui| {
                ui.label(egui::RichText::new(format!("error: {err}")).color(Color32::from_rgb(230, 80, 80)));
                if ui.small_button("dismiss").clicked() {
                    self.error = None;
                }
            });
        }
    }

    fn canvas(&mut self, ui: &mut egui::Ui) {
        let (response, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
        let rect = response.rect;
        painter.rect_filled(rect, 0.0f32, Color32::from_rgb(8, 9, 14));

        // Interaction.
        if response.dragged_by(egui::PointerButton::Primary) {
            self.view.pan(response.drag_delta());
            self.follow = None;
        }
        if response.dragged_by(egui::PointerButton::Secondary) {
            let d = response.drag_delta();
            self.view.yaw += d.x as f64 * 0.01;
            self.view.pitch = (self.view.pitch + d.y as f64 * 0.01).clamp(-1.55, 1.55);
        }
        if response.hovered() {
            let (scroll, zoom) = ui.ctx().input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
            let factor = zoom as f64 * (1.0 + scroll as f64 * 0.002);
            if factor != 1.0 {
                if let Some(mouse) = response.hover_pos() {
                    // Zoom about the cursor: keep the world point under the mouse fixed.
                    let before = Vec2::new(mouse.x - rect.center().x, mouse.y - rect.center().y);
                    self.view.scale *= factor;
                    let after = before * factor as f32;
                    self.view.pan(after - before);
                }
            }
        }
        if let Some(id) = self.follow {
            if let Some(b) = self.sim.system.get(id) {
                self.view.center = b.position;
            } else {
                self.follow = None;
            }
        }

        // Projected bodies, back to front.
        let mut items: Vec<(usize, Pos2, f64)> = self.sim.system.bodies().iter().enumerate().map(|(i, b)| {
            let (p, z) = self.view.project(b.position, rect);
            (i, p, z)
        }).collect();
        items.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        if response.clicked() {
            if let Some(mouse) = response.interact_pointer_pos() {
                let mut best: Option<(f32, usize)> = None;
                for (i, p, _) in &items {
                    let d = p.distance(mouse);
                    if d < 14.0 && best.map(|(bd, _)| d < bd).unwrap_or(true) {
                        best = Some((d, *i));
                    }
                }
                self.selected = best.map(|(_, i)| self.sim.system.bodies()[i].id());
            }
        }

        // Trails.
        if self.show_trails {
            for (id, trail) in &self.trails {
                if trail.len() < 2 {
                    continue;
                }
                let colour = self.sim.system.get(*id).map(|b| b.colour()).unwrap_or([120, 120, 120]);
                let pts: Vec<Pos2> = trail.iter().map(|p| self.view.project(*p, rect).0).collect();
                let c = c32(colour);
                painter.add(egui::Shape::line(pts, Stroke::new(1.0f32, Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), 110))));
            }
        }

        // Bodies.
        let selected = self.selected;
        for (i, p, _) in &items {
            let b = &self.sim.system.bodies()[*i];
            if !rect.expand(50.0).contains(*p) {
                continue;
            }
            let physical = (b.radius() * self.view.scale) as f32;
            let r = physical.max(if b.kind == BodyKind::Star { 5.0 } else { 3.0 });
            let colour = c32(b.colour());
            painter.circle_filled(*p, r, colour);
            if b.kind == BodyKind::Star {
                painter.circle_stroke(*p, r + 3.0, Stroke::new(1.0f32, Color32::from_rgba_unmultiplied(255, 220, 120, 90)));
            }
            if b.kind == BodyKind::Projectile {
                painter.circle_stroke(*p, r + 4.0, Stroke::new(1.5f32, Color32::from_rgb(255, 90, 90)));
            }
            if selected == Some(b.id()) {
                painter.circle_stroke(*p, r + 6.0, Stroke::new(1.5f32, Color32::WHITE));
            }
            if self.show_labels || selected == Some(b.id()) {
                painter.text(*p + Vec2::new(r + 5.0, 0.0), Align2::LEFT_CENTER, &b.name, FontId::proportional(12.0), Color32::from_gray(210));
            }
        }

        // Scale bar.
        let target_px = 120.0;
        let metres = target_px / self.view.scale;
        let unit_len = metres / self.units.length;
        let nice = 10f64.powf(unit_len.log10().floor());
        let nice = [1.0, 2.0, 5.0].iter().map(|m| m * nice).filter(|v| *v <= unit_len).last().unwrap_or(nice);
        let px = (nice * self.units.length * self.view.scale) as f32;
        let p0 = Pos2::new(rect.left() + 16.0, rect.bottom() - 16.0);
        painter.line_segment([p0, p0 + Vec2::new(px, 0.0)], Stroke::new(2.0f32, Color32::GRAY));
        painter.text(p0 + Vec2::new(0.0, -6.0), Align2::LEFT_BOTTOM, format!("{} {}", sig(nice), self.unit_spec.length), FontId::proportional(11.0), Color32::GRAY);
        painter.text(
            Pos2::new(rect.right() - 10.0, rect.bottom() - 10.0),
            Align2::RIGHT_BOTTOM,
            "drag: pan · scroll: zoom · right-drag: tilt · click: select · space: run/pause · F: fit",
            FontId::proportional(11.0),
            Color32::from_gray(90),
        );
        if self.sim.system.is_empty() {
            painter.text(rect.center(), Align2::CENTER_CENTER, "Add bodies on the left or load a preset above", FontId::proportional(16.0), Color32::GRAY);
        }
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let (space, fit, delete) = ctx.input(|i| (i.key_pressed(egui::Key::Space), i.key_pressed(egui::Key::F), i.key_pressed(egui::Key::Delete)));
        if space && !ctx.wants_keyboard_input() {
            self.running = !self.running;
        }
        if fit && !ctx.wants_keyboard_input() {
            self.fit_view();
        }
        if delete && !ctx.wants_keyboard_input() {
            if let Some(id) = self.selected.take() {
                self.sim.system.remove(id);
                self.sim.reset_baseline();
            }
        }
        if self.running {
            let dt = ctx.input(|i| i.stable_dt) as f64;
            self.advance_frame(dt);
            ctx.request_repaint();
        }

        egui::TopBottomPanel::top("top").show(ctx, |ui| self.top_panel(ui));
        egui::TopBottomPanel::bottom("bottom").show(ctx, |ui| self.bottom_panel(ui));
        // Panels take a share of the window so the orbit view stays visible on
        // small screens; on a narrow one the projectile panel starts hidden.
        let width = ctx.screen_rect().width();
        if !self.layout_initialised && width > 0.0 {
            self.layout_initialised = true;
            if width < 1000.0 {
                self.show_projectile = false;
            }
            if width < 640.0 {
                self.show_bodies = false;
            }
        }
        let panel_width = (width * 0.27).clamp(220.0, 360.0);
        if self.show_bodies {
            egui::SidePanel::left("bodies").resizable(true).default_width(panel_width).max_width(width * 0.45).show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| self.bodies_panel(ui));
            });
        }
        if self.show_projectile {
            egui::SidePanel::right("projectile").resizable(true).default_width(panel_width).max_width(width * 0.45).show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| self.projectile_panel(ui));
            });
        }
        egui::CentralPanel::default().show(ctx, |ui| self.canvas(ui));
    }
}

#[cfg(target_arch = "wasm32")]
mod web {
    use super::App;
    use eframe::egui;
    use wasm_bindgen::prelude::*;
    use wasm_bindgen::JsCast;

    /// Start the simulator inside `container`: a `<canvas class="nbody-canvas">`
    /// is created there and the app is drawn into it with WebGL2. The returned
    /// promise resolves once the app is running and rejects (with a JS value)
    /// if the browser lacks WebGL2.
    #[wasm_bindgen]
    pub async fn start(container: web_sys::HtmlElement) -> Result<(), JsValue> {
        console_error_panic_hook::set_once();
        let document = web_sys::window().ok_or("no window")?.document().ok_or("no document")?;
        let canvas: web_sys::HtmlCanvasElement = document.create_element("canvas")?.dyn_into()?;
        canvas.set_class_name("nbody-canvas");
        canvas.set_attribute("aria-label", "N-body simulation: orbit view and controls")?;
        container.append_child(&canvas)?;
        eframe::WebRunner::new()
            .start(
                canvas,
                eframe::WebOptions::default(),
                Box::new(|cc| {
                    cc.egui_ctx.set_visuals(egui::Visuals::dark());
                    Ok(Box::new(App::new()))
                }),
            )
            .await?;
        // eframe adds an off-screen <input> that relays keyboard input to the
        // canvas; give it a name so accessibility checks do not flag it.
        if let Ok(Some(agent)) = document.query_selector("body > input[autocapitalize=\"off\"]:not([aria-label])") {
            agent.set_attribute("aria-label", "Keyboard input for the N-body simulation")?;
        }
        Ok(())
    }
}
