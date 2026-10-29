# n-body-problem: Newtonian gravity for layered, colliding worlds

A Rust simulator for the gravitational N-body problem that lets you build a
planetary system body by body, change each body's mass, radius, material and
internal layering, watch it evolve with symplectic integrators, and fire
projectiles at it from infinity to see whether they are **captured**,
**ejected** or **swallowed**, and whether they knock a planet out on the way.

It ships as three crates:

| crate | what it is |
|-------|------------|
| `nbody-core` | the physics library (no UI dependencies) |
| `nbody` | a command-line driver with an interactive shell |
| `nbody-gui` | an interactive window (`nbody-desktop`) and, compiled to WebAssembly, the same app in a web page |

Everything is Newtonian point-mass gravity with finite-size spheres for
collisions. There is no relativity, no tides, no spin dynamics, no
fragmentation; see [Limitations](#9-limitations) for the honest list.

---

## Contents

1. [Quick start](#1-quick-start)
2. [The physics in plain language](#2-the-physics-in-plain-language)
3. [Projectile mode](#3-projectile-mode)
4. [The numerics](#4-the-numerics)
5. [Scenario files](#5-scenario-files)
6. [Command-line reference](#6-command-line-reference)
7. [The GUI](#7-the-gui)
8. [Tour of the code](#8-tour-of-the-code)
9. [Limitations](#9-limitations)
10. [How we know it works](#10-how-we-know-it-works)

---

## 1. Quick start

You need a Rust toolchain (`rustup`, stable). Then:

```sh
git clone https://github.com/BlankRelax/n-body-problem.git
cd n-body-problem
cargo build --release            # builds nbody (CLI) and nbody-desktop (GUI)

# Integrate the Solar System for 12 years and report energy drift and orbital elements
./target/release/nbody run solar-system

# Fire a rogue Earth-mass planet at Jupiter and see what happens to it
./target/release/nbody projectile projectile-capture

# Map capture / merger / ejection over impact parameter and speed
./target/release/nbody sweep projectile-capture --b 19:24:41 --v-inf 0.5:2:4

# Open the interactive window
./target/release/nbody-desktop
```

`nbody preset` lists the built-in scenarios (Solar System, Earth–Moon, the
figure-eight choreography, Burrau's Pythagorean problem, a binary star with
a planet, a 200-star Plummer cluster, and two projectile experiments).
`nbody repl` opens a shell where you can `add`, `set`, `run`, `projectile`
and `save` interactively.

### In the browser

The GUI also compiles to WebAssembly and runs in any browser with WebGL2
(it is live at [hassaanmirza.co.uk/n-body](https://hassaanmirza.co.uk/n-body)).
You need the `wasm32-unknown-unknown` target and `wasm-bindgen-cli` at
exactly the version pinned in `crates/nbody-gui/Cargo.toml`:

```sh
rustup target add wasm32-unknown-unknown
cargo install wasm-bindgen-cli --version 0.2.129 --locked
node scripts/build.mjs      # dist/ (a static site) and dist-single/n-body-problem.html
```

`dist/` is a static site: serve it from any web server, or embed the module
elsewhere by importing `pkg/nbody_gui.js` and calling its exported
`start(container)`. `dist-single/n-body-problem.html` is one self-contained
file that opens straight from disk. Pages with a Content-Security-Policy
need `'wasm-unsafe-eval'` in `script-src`.

---

## 2. The physics in plain language

**Gravity.** Every body pulls on every other with Newton's inverse-square
force. For `N` bodies that is `N(N−1)/2` pair forces per step, which the
code sums exactly by default and approximates with a Barnes–Hut tree when
`N` gets large (hundreds to hundreds of thousands).

**Bodies are layered spheres.** Each body is a stack of concentric shells,
each made of one material (iron, rock, ice, water, gas, stellar plasma,
carbon, regolith, rubber, or a custom density). A famous result of Newton,
the *shell theorem*, says that the gravity outside any spherically
symmetric body is the same as a point mass at its centre, so the layering
never changes the orbits. What it does fix is

* the **radius** for a given mass, and so when two bodies touch, the speed
  needed to escape from their surfaces, and the Roche limit inside which
  tides would tear a satellite apart;
* the **moment of inertia** (reported as `I/MR²`: 0.4 for a uniform sphere,
  about 0.33 for the Earth, 0.07 for the Sun);
* what a **merger** produces: the materials are re-stacked densest-first,
  the way gravity sorts a molten planet, so iron sinks to the new core and
  ice floats to the surface.

Changing a body's mass rescales it homologously (densities fixed, radius
∝ mass^(1/3)); changing its material at fixed mass recomputes its radius;
you can also add or remove shells one at a time. The Solar System presets
build each planet as a two-layer core + mantle model solved to reproduce its
real mass *and* radius (the fit tells you, for instance, that the Earth
needs a 0.8-Earth-radius iron core if iron and rock are uncompressed).

**Collisions.** When two spheres overlap the simulator either merges them
(mass and momentum conserved exactly), bounces them (an impulsive collision
with the materials' coefficients of restitution), or decides *physically*
using the regime map of Leinhardt & Stewart (2012): slower than the mutual
escape speed they cannot separate again and merge; faster, they merge if the
hit is head-on and bounce off each other (hit-and-run) if it is grazing. The
energy lost in each event is booked separately, so the reported energy drift
is purely numerical.

**Integration.** Orbits are advanced with *symplectic* integrators (leapfrog,
Yoshida's 4th-order scheme and the optimised PEFRL scheme), which conserve
energy to round-off over millions of orbits instead of slowly spiralling
in like ordinary Runge–Kutta (which is included so you can see the
difference). Close encounters can use adaptive steps tied to the shortest
local orbital time.

---

## 3. Projectile mode

An object arriving from far away is described by two numbers: its speed at
infinity `v∞` and its *impact parameter* `b`, the sideways offset of the
straight line it would follow if gravity were switched off. Relative to the
system's centre of mass the incoming path is a hyperbola, and the simulator
starts the projectile *exactly on that hyperbola* at a chosen distance
(default: ten times the size of the system) rather than guessing a launch
velocity. The two-body prediction for the closest approach
`q = sqrt(a² + b²) − a` with `a = GM/v∞²` and the Rutherford deflection
angle `2 atan(GM/(b v∞²))` are printed for comparison.

A single star cannot capture anything: energy conservation sends the
visitor back out with the same `v∞`. Capture needs a planet to absorb
energy, and the same exchange can kick the planet out instead. The tracker
therefore watches the projectile's energy relative to the rest of the system
and reports

| verdict | meaning |
|---------|---------|
| **captured** | bound (ε < 0) after clearing every planet's Hill sphere; the osculating orbit is printed |
| **ejected** | unbound and receding past its start distance; the outgoing `v∞` and deflection are printed |
| **merged** | absorbed by a body (the collision log says which, at what speed and angle) |
| plus | any original body that the visit left unbound, and the closest surface-to-surface miss |

The built-in `projectile-capture` scenario aims a rogue Earth-mass planet so
that its periapsis lies on Jupiter's orbit and phases Jupiter to be a chosen
angle from that point at the analytically computed arrival time. Scanning
that angle reproduces the textbook regimes of a gravity assist (output of
`cargo test -p nbody-core scan_capture_lag -- --ignored --nocapture`,
abridged):

```
lag -4.0°: captured   closest gap 0.048 au   bound orbit a = 18 au, e = 0.73
lag -2.5°: captured   closest gap 0.012 au   bound orbit a = 5.6 au, e = 0.44   <- default preset
lag -1.2°: captured   closest gap 3e-4  au   bound orbit a = 3.1 au, e = 0.69
lag  0.0°: merged     hits Jupiter
lag +1.2°: captured   closest gap 4e-4  au   bound orbit a = 3.3 au, e = 0.64
lag +3.6°: ejected    leaves at 3.3 km/s (arrived at 4.7 km/s: braked, but not enough)
lag +4.5°: ejected    leaves at 5.9 km/s (slingshot: gained energy)
lag +6.0°: ejected    leaves at 6.6 km/s
```

`nbody sweep` does the same over a grid of `(v∞, b)` and prints the map,
with each outcome's share of the cross-section (weighted by `2πb db`).

---

## 4. The numerics

* **Force:** softened pair force `−G m_j (r_i − r_j)/(|r_ij|² + ε²)^{3/2}`
  (ε = 0 by default; useful for collisionless clusters). Direct summation or
  Barnes–Hut octree with opening angle θ (θ = 0 reproduces direct summation
  exactly; θ = 0.5 gives ~1 % force errors).
* **Integrators:** leapfrog KDK (2nd order, 1 force evaluation per step),
  Yoshida 4 (3 evaluations), PEFRL (4 evaluations, smallest error constant),
  RK4 (4 evaluations, non-symplectic). Fixed step or adaptive
  `Δt = η · min_ij sqrt(r_ij³ / G(m_i + m_j))`.
* **Conservation tracking:** energy, momentum and angular momentum drift
  relative to a baseline, net of collision bookkeeping.
* **Units:** SI inside; scenario files pick any unit system (au / solar
  masses / years, km / kg / hours, pc / solar masses / Myr, ...). `G` can be
  overridden for textbook problems with `G = 1`.

Measured on the presets (release build, defaults):

| scenario | integrator | step | run | max |ΔE/E| |
|----------|-----------|------|-----|------------|
| Solar System, 9 bodies | Yoshida 4 | 0.0005 yr | 12 yr | 2e-11 |
| figure-eight, G = 1 | PEFRL | 0.002 | 3 periods | 1.5e-12 |
| Pythagorean problem (close encounters) | PEFRL, adaptive | ≤ 0.01 | t = 70 | 1e-8 |
| rogue planet past Jupiter | PEFRL, adaptive | ≤ 0.002 yr | 60 yr | 3e-13 |
| Plummer cluster, 200 stars | leapfrog + Barnes–Hut θ = 0.5 | 0.002 Myr | 1 Myr | 1.3e-4 |

Details and references are in [docs/physics.md](docs/physics.md).

---

## 5. Scenario files

A scenario is a TOML file. Angles are in degrees, custom densities always in
kg/m³, everything else in the declared units. Bodies are positioned either
explicitly or by Keplerian elements relative to a body listed above them.

```toml
name = "Sun, Earth and a comet"

[units]
length = "au"          # m, km, au, ly, pc, r_sun, r_earth
mass = "solar"         # kg, sun/solar, earth, jupiter, moon
time = "year"          # s, min, hour, day, year, myr

[simulation]
integrator = "yoshida4"     # leapfrog | yoshida4 | pefrl | rk4
dt = 0.001                  # or: [simulation.adaptive] eta = 0.01, dt_min = 1e-8, dt_max = 0.01
collisions = "physical"     # none | merge | bounce | physical
force = "direct"            # direct | barnes_hut (with theta = 0.5)
softening = 0.0
duration = 10.0
# g = 1.0                   # override G in these units (textbook problems)

[[bodies]]
name = "Sun"
kind = "star"
mass = 1.0
material = "plasma"         # mass + material -> radius from density

[[bodies]]
name = "Earth"
mass = 3.003e-6
radius = 4.26e-5
material = "rock"           # mass + radius + material + core -> two-layer fit
core = "iron"
orbit = { around = "Sun", a = 1.0, e = 0.0167, inclination = 0.0, mean_anomaly = 90.0 }

[[bodies]]
name = "Layered moon"
position = [1.0026, 0, 0]
velocity = [0, 6.5, 0]
layers = [                  # explicit shells, inside out
  { outer_radius = 5e-6, material = "iron" },
  { outer_radius = 1.1e-5, material = "rock" },
  { outer_radius = 1.2e-5, material = "ice" },
]

[projectile]
name = "Comet"
radius = 3.3e-8
material = "ice"
v_infinity = 3.0            # au/year
impact_parameter = 0.3      # au
direction = [1, 0, 0]       # velocity direction at infinity
max_time = 30.0
```

`nbody preset <name> -o file.toml` writes any preset out as a starting point;
the `scenarios/` directory holds all of them. `nbody run ... --save` and the
shell's `save` command snapshot a live system with explicit states and
layers.

---

## 6. Command-line reference

```
nbody run <scenario|preset> [--duration T] [--dt] [--integrator] [--collisions]
                            [--output traj.csv] [--samples N] [--save state.toml]
nbody projectile <scenario|preset> [--v-inf] [--b] [--mass] [--radius] [--material]
                            [--core] [--direction x,y,z] [--offset x,y,z]
                            [--distance] [--max-time] [--output traj.csv]
nbody sweep <scenario|preset> --b lo:hi:n --v-inf lo:hi:n [--log-b] [--max-time] [--json out.json]
nbody info <scenario|preset>       bodies, osculating elements, Hill radii, Roche limits, suggested dt
nbody preset [name] [-o file]      list presets or write one out
nbody materials                    the material table
nbody repl [scenario|preset]       interactive shell
```

The shell's commands (type `help`):

```
add Ceres mass=4.7e-10 material=rock orbit=Sun,2.77,0.08,10.6
set Mars mass=6.4e-7 material=ice
shell Earth add 2e-5 ice            # add an icy shell
circular Moon around Earth
orbit Probe around Earth a=6e-5 e=0.3 i=28
run 5                               # 5 time units
projectile vinf=2 b=0.4 mass=1e-9 material=ice dir=1,0,0
outcome
sweep vinf=0.5:2:4 b=0.1:3:30
save my_system.toml
```

---

## 7. The GUI

`nbody-desktop` (or the web build at `/n-body`) opens with the Solar System loaded.

* **Top bar:** choose a preset or open/save a scenario file; run / pause /
  step / reset; simulation speed (log slider); integrator; collision mode;
  time step or adaptive stepping; barycentre; fit view; labels and trails.
* **Left panel:** the body list. Select one to edit its name, mass, radius,
  material, kind, position and velocity, inspect and edit its shells, put it
  on a circular orbit around another body, follow it with the camera, or
  delete it. The "Add a body" section creates a new body either on a
  Keplerian orbit (semi-major axis, eccentricity, inclination) about any
  existing body or at an explicit position and velocity.
* **Right panel:** projectile mode. Set `v∞`, `b`, mass, material, approach
  direction and which side of the system to pass on; the two-body periapsis
  and deflection are predicted live. "Launch" injects it and runs; the
  verdict, closest approach, closest miss distance, energy history and any
  displaced bodies update as it flies. The outcome sweep runs a grid of
  experiments on a copy of the system.
* **Bottom:** clock, step counts, energy and virial ratio, relative drifts
  of energy / momentum / angular momentum, an energy-error sparkline and the
  event log (collisions, launches).
* **Canvas:** drag to pan, scroll to zoom, right-drag to tilt into 3-D,
  click to select, space to run/pause, F to fit, Delete to remove.

---

## 8. Tour of the code

```
crates/nbody-core/src
  vector.rs        3-vector
  units.rs         constants, unit systems, formatting
  material.rs      densities, restitution, colours
  structure.rs     layered spheres: mass, radius, inertia, enclosed mass, two-layer fit, differentiated merger
  body.rs          a body and its edits (set_mass / set_radius / set_material)
  system.rs        the body collection, G, softening, clock; circular / Keplerian placement
  force/           direct summation and Barnes-Hut octree; softened potential; dynamical time
  integrator/      leapfrog, Yoshida 4, PEFRL, RK4
  collision.rs     merge / bounce / physical regimes with energy bookkeeping
  orbit.rs         elements <-> state, Kepler's equation, Hill / Roche / SOI, hyperbolic flybys
  diagnostics.rs   energy, momentum, angular momentum, virial ratio; conservation tracker
  projectile.rs    hyperbolic injection, approach geometry, outcome tracker, sweeps
  simulation.rs    the driver: time-step policy, stepping, collisions, trackers
  scenario.rs      TOML format, build to a Simulation, snapshots
  presets.rs       built-in scenarios (incl. the phased Jupiter encounter)
crates/nbody-cli/src
  main.rs          clap subcommands     repl.rs  interactive shell     report.rs  tables
crates/nbody-gui/src
  main.rs          eframe/egui application
```

The core is `#![forbid(unsafe_code)]` and depends only on `serde`, `toml`,
`rand` and `thiserror`.

---

## 9. Limitations

* Point-mass Newtonian gravity only: no post-Newtonian terms, no tides, no
  oblateness, no radiation pressure or drag.
* Bodies do not rotate; the moment of inertia is reported but unused by the
  dynamics.
* Collisions never fragment. A fast head-on impact accretes where reality
  would shatter both bodies.
* Mergers assume complete differentiation (densest material innermost).
* The Barnes–Hut tree uses monopoles only and is rebuilt every evaluation;
  it is single-threaded.
* Adaptive stepping is global (one Δt for all bodies) and breaks exact
  symplecticity; use it for few-body encounters, not for long cluster runs.
* The GUI runs the integration on the UI thread and caps the work per frame,
  so very expensive settings will slow the display rather than the clock.

---

## 10. How we know it works

`cargo test` runs 61 unit tests, among them:

* **Force = −∇U** by finite differences on a softened 5-body configuration.
* **Barnes–Hut vs direct summation:** identical at θ = 0, RMS error < 1 % at
  θ = 0.5 on a 2000-body cloud; coincident bodies terminate.
* **Integrator orders** measured on an eccentric Kepler orbit (≥ 2 for
  leapfrog, ≥ 4 for the others); symplectic schemes keep bounded energy error
  over 200 orbits while RK4 drifts; PEFRL beats Yoshida at equal step.
* **Earth returns to its starting point** after one year with position error
  < 1e-6 au and energy error < 1e-10.
* **Kepler element round trips** (elliptic and hyperbolic, all angles) and
  both forms of Kepler's equation to 1e-12.
* **Solar System preset periods** match `a^{3/2}` for all eight planets.
* **Hyperbolic injection** has exactly the requested energy and angular
  momentum; the analytic time and point of periapsis agree with a numerical
  integration; a single-body flyby is classified *ejected* with the Rutherford
  deflection angle to 0.2 %; a head-on shot *merges* with exact mass
  conservation.
* **Collision laws:** merger conserves mass and momentum and changes the
  energy by exactly `−½μv² + G m₁m₂/d`; bounce conserves momentum with the
  prescribed restitution; the physical regime map merges slow impacts and
  bounces fast grazing ones.
* **Jupiter encounter regimes:** the default lag captures, zero lag merges,
  +4.5° ejects.
* The figure-eight choreography returns to its initial state after one
  period; the Plummer sampler produces a cluster in virial equilibrium.
