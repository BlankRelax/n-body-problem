# Physics and numerics notes

This document records the equations the code implements and the choices
behind them, with pointers to the source. It assumes a reader comfortable
with classical mechanics; the README has the plain-language version.

## 1. Equations of motion

For `N` point masses the Hamiltonian is

```
H = Σ_i p_i² / (2 m_i)  −  Σ_{i<j} G m_i m_j / |q_i − q_j|
```

so that

```
a_i = −G Σ_{j≠i} m_j (q_i − q_j) / |q_i − q_j|³
```

With a Plummer softening length `ε` the pair term becomes
`−G m_i m_j / sqrt(r² + ε²)`, whose gradient is the force actually evaluated
(`force/mod.rs`, `force/direct.rs`). The force and potential are therefore
exactly conjugate and `K + U` is an exact invariant of the softened system,
which is what `diagnostics.rs` monitors. `ε = 0` is the default for planetary
systems: finite-size collisions are handled geometrically instead, so the
singularity is never reached.

**Shell theorem.** The gravitational field outside a spherically symmetric
mass distribution equals that of a point mass at its centre. Every body here
is a stack of concentric uniform shells (`structure.rs`), so the layered
structure changes only its radius, moment of inertia and merger behaviour,
never the orbital dynamics. This is exact, not an approximation.

## 2. Integrators (`integrator/`)

The Hamiltonian is separable, `H = T(p) + U(q)`, so the flow can be split
into exact kicks `p ← p − ∇U Δt` and drifts `q ← q + p/m Δt`. Compositions of
these are *symplectic maps*: they preserve phase-space volume and, by
backward-error analysis, exactly conserve a modified Hamiltonian
`H̃ = H + O(Δt^k)`. The energy error therefore stays bounded for
exponentially long times rather than growing secularly, which is why
symplectic schemes are the standard for long orbital integrations
(Wisdom & Holman 1991; Hairer, Lubich & Wanner 2006).

| scheme | order | force evals / step | notes |
|--------|-------|--------------------|-------|
| Leapfrog KDK | 2 | 1 | time-reversible; the workhorse for collisionless N-body |
| Yoshida (1990) | 4 | 3 | triple-jump composition `S₂(w₁h) S₂(w₀h) S₂(w₁h)`, `w₁ = 1/(2−2^{1/3})`, `w₀ = 1 − 2w₁`; the middle stage runs backwards in time |
| PEFRL (Omelyan, Mryglod & Folk 2002) | 4 | 4 | position-extended Forest–Ruth-like; coefficients optimised for the smallest error constant |
| RK4 | 4 | 4 | not symplectic; included as the control: its energy error grows linearly |

`integrator/mod.rs` contains convergence tests that verify the observed
orders on an eccentric Kepler orbit and that RK4 drifts while the symplectic
schemes do not.

**Adaptive steps.** `TimestepPolicy::Adaptive` sets
`Δt = η · min_{i<j} sqrt((r_ij² + ε²)^{3/2} / (G (m_i + m_j)))`, i.e. a
fraction of the shortest local two-body dynamical time `1/ω`. A state-
dependent step breaks exact symplecticity (the map is no longer the flow of a
fixed modified Hamiltonian), but for the few-body close encounters of
projectile mode it is far better than a fixed step that is either wasteful
far away or unresolved at periapsis. The tests confirm the energy error of
such runs stays at the 1e-8 level over a 120-year rogue-planet encounter.

## 3. Barnes–Hut tree (`force/barnes_hut.rs`)

An octree is rebuilt at every force evaluation. Each node stores its total
mass and centre of mass; a body interacts with a node as a single monopole
when the opening criterion `s / d < θ` holds (`s` the node width, `d` the
distance to its centre of mass), otherwise the node's children are visited.
Cost is `O(N log N)` with a controllable error: with `θ = 0.5` the RMS
relative force error on a random cloud is ~1e-2, with `θ = 0` the result is
identical to direct summation (tested). Because the tree force is not exactly
antisymmetric, momentum is conserved only to the force-error level; the
Plummer-cluster preset shows this as a ~1e-4 momentum drift per crossing
time.

## 4. Collisions (`collision.rs`)

Two spheres collide when `|r_ij| < R_i + R_j`. Three resolutions exist:

* **Merge** conserves mass and momentum exactly. The remnant is placed at the
  centre of mass; its structure is the *differentiated* union of both bodies:
  every material keeps its volume and the shells are re-stacked densest-first,
  so iron sinks to the core and ice floats to the surface. The system's
  mechanical energy `K + U` changes by the lost relative kinetic energy
  `−½ μ v_rel²` (μ the reduced mass) *plus* the removed pair potential
  `+G m_i m_j / d`, which becomes internal binding energy of the remnant.
* **Bounce** applies the impulse `J = −(1 + e) (v_rel·n) / (1/m_i + 1/m_j)`
  along the contact normal `n`, with `e = sqrt(e_i e_j)` the geometric mean of
  the materials' coefficients of restitution, then separates the spheres.
* **Physical** (default) follows the regime map of Leinhardt & Stewart (2012):
  impacts slower than the mutual escape speed `sqrt(2G(m_i+m_j)/(R_i+R_j))`
  always merge; faster impacts merge when non-grazing and hit-and-run (bounce)
  when grazing, with the grazing boundary `sin θ > R_target / (R_target +
  R_projectile)`. Catastrophic disruption and fragmentation are not modelled.

The conservation tracker books every event's `ΔE` so the reported energy
error is purely numerical.

## 5. Keplerian machinery (`orbit.rs`)

State ↔ element conversions use the standard perifocal frame with the
rotation `R_z(Ω) R_x(i) R_z(ω)`. Circular and equatorial degeneracies are
handled explicitly (argument of latitude / true longitude). Kepler's
equation is solved by Newton iteration in both the elliptic (`E − e sin E =
M`) and hyperbolic (`e sinh H − H = M`) forms. Scale lengths provided:

* Hill radius `a(1−e)(m/3M)^{1/3}`,
* Laplace sphere of influence `a (m/M)^{2/5}`,
* rigid Roche limit `R_p (2ρ_p/ρ_s)^{1/3}` and fluid limit `2.44 R_p (ρ_p/ρ_s)^{1/3}`,
* hyperbolic periapsis `q = sqrt(a² + b²) − a` with `a = μ/v∞²`,
* Rutherford deflection `θ = 2 atan(μ / (b v∞²))`,
* gravitationally focused cross-section `π R² (1 + v_esc²/v∞²)`.

## 6. Projectile mode (`projectile.rs`)

An interloper is specified by its excess speed `v∞`, impact parameter `b`
and approach direction. Relative to the system's barycentre and total mass
the incoming leg is a two-body hyperbola with

```
a = μ/v∞²,   e = sqrt(1 + b²/a²),   p = b²/a,   q = p/(1+e)
```

`hyperbolic_approach_state` evaluates it at the chosen start distance `d`:
the true anomaly from `r = p/(1 + e cos ν)` on the incoming branch, the
velocity from the vis-viva relations `v_r = sqrt(μ/p) e sin ν`,
`v_⊥ = sqrt(μ/p)(1 + e cos ν)`, and a rotation that maps the asymptote
direction onto the requested approach direction and the asymptote's offset
`r_⊥ = v̂ × h / |v|` onto the requested offset direction. The hyperbolic
Kepler equation then gives the flight time to periapsis
(`approach_geometry`), which `presets::projectile_capture_with_lag` uses to
phase Jupiter so that it is at a chosen angle from the periapsis point when
the projectile arrives.

**Why capture needs a third body.** In the two-body problem energy is
conserved along the hyperbola, so a body arriving with `v∞ > 0` leaves with
the same `v∞`. Only an exchange of energy with another body (a planet) can
drive the specific energy

```
ε = ½ |v − v_cm|² − Σ_j G m_j / |r − r_j|
```

negative. The tracker therefore

* suspends all verdicts while the projectile is within three Hill radii of
  any secondary (its energy is then dominated by that body's local well and
  its barycentric radial motion is meaningless);
* declares **Captured** once it is past periapsis, clear of every secondary
  and `ε < 0`: energy changes only in encounters, so it will come back;
* declares **Ejected** once `ε > 0`, it is receding and it has passed its
  start distance, reporting the outgoing `v∞` and the deflection angle;
* reports **Merged** from the collision log;
* lists which original members were themselves unbound by the visit.

The scan in `presets.rs` (`scan_capture_lag`) shows the three regimes as
Jupiter's phase is varied: a pass in front of Jupiter brakes the rogue
planet into a bound heliocentric orbit, a direct hit merges, and a pass
behind slingshots it out faster than it arrived. The regression test
`jupiter_encounter_regimes` pins one representative of each.

## 7. Units

Everything internal is SI. Scenario files declare a unit system (for
instance au, solar masses and years, in which `G ≈ 4π²`) and are converted
on load; outputs are converted back. The gravitational constant is a field
of the system, so textbook problems with `G = 1` are expressed by overriding
it in the scenario.

## References

* Barnes, J. & Hut, P. (1986), *Nature* 324, 446.
* Aarseth, S. J., Hénon, M. & Wielen, R. (1974), *A&A* 37, 183 (Plummer sampling).
* Chenciner, A. & Montgomery, R. (2000), *Ann. Math.* 152, 881 (figure eight).
* Hairer, E., Lubich, C. & Wanner, G. (2006), *Geometric Numerical Integration*, Springer.
* Leinhardt, Z. M. & Stewart, S. T. (2012), *ApJ* 745, 79 (collision regimes).
* Omelyan, I. P., Mryglod, I. M. & Folk, R. (2002), *Comput. Phys. Commun.* 146, 188 (PEFRL).
* Szebehely, V. & Peters, C. F. (1967), *AJ* 72, 876 (Pythagorean problem).
* Wisdom, J. & Holman, M. (1991), *AJ* 102, 1528 (symplectic maps for N-body).
* Yoshida, H. (1990), *Phys. Lett. A* 150, 262.
