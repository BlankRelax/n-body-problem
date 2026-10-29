//! Two-body (Keplerian) orbital mechanics.
//!
//! Everything here is the analytic solution of the two-body problem with
//! gravitational parameter `μ = G(M + m)`. It is used to *place* bodies
//! (orbital elements → state vector), to *describe* their motion (state
//! vector → osculating elements), and for the classical scale lengths of
//! planetary dynamics (Hill sphere, Roche limit, sphere of influence) that
//! tell you whether a configuration is stable before you simulate it.

use crate::vector::Vec3;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum OrbitError {
    #[error("gravitational parameter must be positive (got {0})")]
    NonPositiveMu(f64),
    #[error("eccentricity must be non-negative (got {0})")]
    BadEccentricity(f64),
    #[error("semi-major axis must be non-zero and finite (got {0})")]
    BadSemiMajorAxis(f64),
    #[error("exactly parabolic orbits (e = 1) have no semi-major axis; use e = 1 ± 1e-9")]
    Parabolic,
    #[error("true anomaly {0} rad is beyond the asymptote of this hyperbola")]
    BeyondAsymptote(f64),
}

/// Classical Keplerian elements, SI and radians.
///
/// For hyperbolic orbits (`e > 1`) the semi-major axis is negative by the
/// usual convention; a positive value is accepted and interpreted as `|a|`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Elements {
    /// Semi-major axis, m.
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination to the reference (x-y) plane, rad.
    pub inclination: f64,
    /// Longitude of the ascending node Ω, rad.
    pub raan: f64,
    /// Argument of periapsis ω, rad.
    pub arg_periapsis: f64,
    /// True anomaly ν, rad.
    pub true_anomaly: f64,
}

impl Elements {
    pub fn circular(a: f64) -> Elements {
        Elements { a, e: 0.0, inclination: 0.0, raan: 0.0, arg_periapsis: 0.0, true_anomaly: 0.0 }
    }

    /// Semi-latus rectum `p = a(1 − e²)` (always positive).
    pub fn semi_latus_rectum(&self) -> f64 {
        (self.a * (1.0 - self.e * self.e)).abs()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OrbitKind {
    Elliptic,
    Parabolic,
    Hyperbolic,
}

/// Osculating elements plus the derived quantities people actually ask for.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OrbitalElements {
    pub elements: Elements,
    pub kind: OrbitKind,
    /// Specific orbital energy `v²/2 − μ/r`, J/kg.
    pub specific_energy: f64,
    /// Specific angular momentum vector `r × v`.
    pub h: Vec3,
    /// Periapsis distance, m.
    pub periapsis: f64,
    /// Apoapsis distance, m (bound orbits only).
    pub apoapsis: Option<f64>,
    /// Orbital period, s (bound orbits only).
    pub period: Option<f64>,
    /// Hyperbolic excess speed, m/s (unbound orbits only).
    pub v_infinity: Option<f64>,
}

impl OrbitalElements {
    pub fn is_bound(&self) -> bool {
        self.specific_energy < 0.0
    }
}

/// Rotation from the perifocal frame to the reference frame: `R_z(Ω) R_x(i) R_z(ω)`.
fn perifocal_to_inertial(el: &Elements, v: Vec3) -> Vec3 {
    let (sw, cw) = el.arg_periapsis.sin_cos();
    let (si, ci) = el.inclination.sin_cos();
    let (so, co) = el.raan.sin_cos();
    // R_z(ω)
    let v = Vec3::new(cw * v.x - sw * v.y, sw * v.x + cw * v.y, v.z);
    // R_x(i)
    let v = Vec3::new(v.x, ci * v.y - si * v.z, si * v.y + ci * v.z);
    // R_z(Ω)
    Vec3::new(co * v.x - so * v.y, so * v.x + co * v.y, v.z)
}

/// Position and velocity relative to the focus for the given elements.
pub fn state_from_elements(mu: f64, el: &Elements) -> Result<(Vec3, Vec3), OrbitError> {
    if !(mu > 0.0) {
        return Err(OrbitError::NonPositiveMu(mu));
    }
    if !(el.e >= 0.0) {
        return Err(OrbitError::BadEccentricity(el.e));
    }
    if el.a == 0.0 || !el.a.is_finite() {
        return Err(OrbitError::BadSemiMajorAxis(el.a));
    }
    if (el.e - 1.0).abs() < 1e-12 {
        return Err(OrbitError::Parabolic);
    }
    let p = el.semi_latus_rectum();
    let nu = el.true_anomaly;
    let denom = 1.0 + el.e * nu.cos();
    if denom <= 0.0 {
        return Err(OrbitError::BeyondAsymptote(nu));
    }
    let r = p / denom;
    let r_pf = Vec3::new(r * nu.cos(), r * nu.sin(), 0.0);
    let k = (mu / p).sqrt();
    let v_pf = Vec3::new(-k * nu.sin(), k * (el.e + nu.cos()), 0.0);
    Ok((perifocal_to_inertial(el, r_pf), perifocal_to_inertial(el, v_pf)))
}

/// Osculating elements of a body with relative position `r` and velocity `v`.
pub fn elements_from_state(mu: f64, r: Vec3, v: Vec3) -> OrbitalElements {
    let tiny = 1e-11;
    let rn = r.norm();
    let v2 = v.norm_sq();
    let h = r.cross(v);
    let hn = h.norm();
    let n = Vec3::Z.cross(h);
    let nn = n.norm();
    let e_vec = (r * (v2 - mu / rn) - v * r.dot(v)) / mu;
    let e = e_vec.norm();
    let energy = 0.5 * v2 - mu / rn;
    let p = hn * hn / mu;

    let kind = if (e - 1.0).abs() < tiny {
        OrbitKind::Parabolic
    } else if e < 1.0 {
        OrbitKind::Elliptic
    } else {
        OrbitKind::Hyperbolic
    };
    let a = match kind {
        OrbitKind::Parabolic => f64::INFINITY,
        _ => -mu / (2.0 * energy),
    };

    let inclination = if hn > 0.0 { (h.z / hn).clamp(-1.0, 1.0).acos() } else { 0.0 };
    let equatorial = nn < tiny * hn.max(1.0);
    let circular = e < tiny;

    let raan = if equatorial { 0.0 } else { n.y.atan2(n.x).rem_euclid(2.0 * PI) };

    let arg_periapsis = if circular {
        0.0
    } else if equatorial {
        // Longitude of periapsis measured from +x; sign from the orbit's sense.
        let w = e_vec.y.atan2(e_vec.x);
        (if h.z >= 0.0 { w } else { -w }).rem_euclid(2.0 * PI)
    } else {
        let mut w = (n.dot(e_vec) / (nn * e)).clamp(-1.0, 1.0).acos();
        if e_vec.z < 0.0 {
            w = 2.0 * PI - w;
        }
        w
    };

    let true_anomaly = if circular {
        // Argument of latitude (from the node) or true longitude (equatorial).
        if equatorial {
            let l = r.y.atan2(r.x);
            (if h.z >= 0.0 { l } else { -l }).rem_euclid(2.0 * PI)
        } else {
            let mut u = (n.dot(r) / (nn * rn)).clamp(-1.0, 1.0).acos();
            if r.z < 0.0 {
                u = 2.0 * PI - u;
            }
            u
        }
    } else {
        let mut nu = (e_vec.dot(r) / (e * rn)).clamp(-1.0, 1.0).acos();
        if r.dot(v) < 0.0 {
            nu = 2.0 * PI - nu;
        }
        nu
    };

    let periapsis = p / (1.0 + e);
    let (apoapsis, period, v_infinity) = match kind {
        OrbitKind::Elliptic => (Some(p / (1.0 - e)), Some(2.0 * PI * (a * a * a / mu).sqrt()), None),
        OrbitKind::Parabolic => (None, None, Some(0.0)),
        OrbitKind::Hyperbolic => (None, None, Some((2.0 * energy).sqrt())),
    };

    OrbitalElements {
        elements: Elements { a, e, inclination, raan, arg_periapsis, true_anomaly },
        kind,
        specific_energy: energy,
        h,
        periapsis,
        apoapsis,
        period,
        v_infinity,
    }
}

// ---- Kepler's equation -------------------------------------------------------

/// Solve `E − e sin E = M` for the eccentric anomaly (elliptic case) by
/// Newton–Raphson with a starting guess that converges for every `e < 1`.
pub fn solve_kepler_elliptic(mean_anomaly: f64, e: f64) -> f64 {
    let m = mean_anomaly.rem_euclid(2.0 * PI);
    let mut ecc = if e < 0.8 { m } else { PI };
    for _ in 0..60 {
        let f = ecc - e * ecc.sin() - m;
        let fp = 1.0 - e * ecc.cos();
        let d = f / fp;
        ecc -= d;
        if d.abs() < 1e-15 {
            break;
        }
    }
    ecc
}

/// Solve `e sinh H − H = M` for the hyperbolic anomaly.
pub fn solve_kepler_hyperbolic(mean_anomaly: f64, e: f64) -> f64 {
    let m = mean_anomaly;
    let mut h = if m.abs() > 1.0 { (m / e).asinh() + (m.abs() / e).ln().max(0.0) * m.signum() } else { m / (e - 1.0).max(1e-3) };
    if !h.is_finite() {
        h = m.signum();
    }
    for _ in 0..100 {
        let f = e * h.sinh() - h - m;
        let fp = e * h.cosh() - 1.0;
        let d = f / fp;
        h -= d;
        if d.abs() < 1e-14 * h.abs().max(1.0) {
            break;
        }
    }
    h
}

/// True anomaly from eccentric anomaly (elliptic).
pub fn true_from_eccentric(ecc_anomaly: f64, e: f64) -> f64 {
    2.0 * ((1.0 + e).sqrt() * (0.5 * ecc_anomaly).sin()).atan2((1.0 - e).sqrt() * (0.5 * ecc_anomaly).cos())
}

/// True anomaly from hyperbolic anomaly.
pub fn true_from_hyperbolic(h: f64, e: f64) -> f64 {
    2.0 * ((e + 1.0).sqrt() * (0.5 * h).tanh()).atan2((e - 1.0).sqrt())
}

/// Mean anomaly → true anomaly for any non-parabolic eccentricity.
pub fn mean_to_true_anomaly(mean_anomaly: f64, e: f64) -> f64 {
    if e < 1.0 {
        true_from_eccentric(solve_kepler_elliptic(mean_anomaly, e), e)
    } else {
        true_from_hyperbolic(solve_kepler_hyperbolic(mean_anomaly, e), e)
    }
}

// ---- scale lengths and speeds ----------------------------------------------------

pub fn circular_velocity(mu: f64, r: f64) -> f64 {
    (mu / r).sqrt()
}

pub fn escape_velocity(mu: f64, r: f64) -> f64 {
    (2.0 * mu / r).sqrt()
}

/// Kepler's third law: `T = 2π sqrt(a³/μ)`.
pub fn kepler_period(mu: f64, a: f64) -> f64 {
    2.0 * PI * (a.abs().powi(3) / mu).sqrt()
}

/// Semi-major axis for a given period.
pub fn semi_major_axis_for_period(mu: f64, period: f64) -> f64 {
    (mu * (period / (2.0 * PI)).powi(2)).cbrt()
}

/// Hill radius of a secondary of mass `m` on an orbit `(a, e)` about `big_m`,
/// evaluated at periapsis: `a(1 − e) (m / 3M)^{1/3}`.
pub fn hill_radius(a: f64, e: f64, m: f64, big_m: f64) -> f64 {
    a * (1.0 - e) * (m / (3.0 * big_m)).cbrt()
}

/// Laplace sphere of influence `a (m/M)^{2/5}`.
pub fn sphere_of_influence(a: f64, m: f64, big_m: f64) -> f64 {
    a * (m / big_m).powf(0.4)
}

/// Rigid-body Roche limit: a satellite of density `rho_s` is torn apart
/// inside `R_p (2 ρ_p / ρ_s)^{1/3}` of a primary of radius `R_p`, density `ρ_p`.
pub fn roche_limit_rigid(primary_radius: f64, rho_primary: f64, rho_satellite: f64) -> f64 {
    primary_radius * (2.0 * rho_primary / rho_satellite).cbrt()
}

/// Fluid Roche limit `≈ 2.44 R_p (ρ_p / ρ_s)^{1/3}`.
pub fn roche_limit_fluid(primary_radius: f64, rho_primary: f64, rho_satellite: f64) -> f64 {
    2.44 * primary_radius * (rho_primary / rho_satellite).cbrt()
}

// ---- hyperbolic encounters --------------------------------------------------------

/// Periapsis distance of a hyperbolic flyby with excess speed `v_inf` and
/// impact parameter `b`: `q = sqrt(a² + b²) − a` with `a = μ / v_inf²`.
pub fn hyperbolic_periapsis(mu: f64, v_inf: f64, b: f64) -> f64 {
    let a = mu / (v_inf * v_inf);
    (a * a + b * b).sqrt() - a
}

/// Impact parameter that produces periapsis `q`: `b = q sqrt(1 + 2μ/(q v_inf²))`
/// (gravitational focusing).
pub fn impact_parameter_for_periapsis(mu: f64, v_inf: f64, q: f64) -> f64 {
    q * (1.0 + 2.0 * mu / (q * v_inf * v_inf)).sqrt()
}

/// Rutherford deflection angle of a two-body hyperbolic flyby:
/// `θ = 2 atan(μ / (b v_inf²))`.
pub fn hyperbolic_deflection(mu: f64, v_inf: f64, b: f64) -> f64 {
    2.0 * (mu / (b * v_inf * v_inf)).atan()
}

/// Gravitationally focused collision cross-section `π R² (1 + v_esc²/v_inf²)`.
pub fn focused_cross_section(mu: f64, radius: f64, v_inf: f64) -> f64 {
    PI * radius * radius * (1.0 + 2.0 * mu / (radius * v_inf * v_inf))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn circular_orbit_round_trip() {
        let mu = 3.986e14;
        let el = Elements::circular(7.0e6);
        let (r, v) = state_from_elements(mu, &el).unwrap();
        assert_relative_eq!(r.norm(), 7.0e6, max_relative = 1e-12);
        assert_relative_eq!(v.norm(), circular_velocity(mu, 7.0e6), max_relative = 1e-12);
        let back = elements_from_state(mu, r, v);
        assert_relative_eq!(back.elements.a, 7.0e6, max_relative = 1e-9);
        assert!(back.elements.e < 1e-9);
        assert_eq!(back.kind, OrbitKind::Elliptic);
        assert_relative_eq!(back.period.unwrap(), kepler_period(mu, 7.0e6), max_relative = 1e-9);
    }

    #[test]
    fn general_elliptic_round_trip() {
        let mu = 1.327e20;
        let el = Elements {
            a: 2.3e11,
            e: 0.37,
            inclination: 0.42,
            raan: 1.1,
            arg_periapsis: 2.9,
            true_anomaly: 4.0,
        };
        let (r, v) = state_from_elements(mu, &el).unwrap();
        let back = elements_from_state(mu, r, v).elements;
        assert_relative_eq!(back.a, el.a, max_relative = 1e-9);
        assert_relative_eq!(back.e, el.e, max_relative = 1e-9);
        assert_relative_eq!(back.inclination, el.inclination, max_relative = 1e-9);
        assert_relative_eq!(back.raan, el.raan, max_relative = 1e-9);
        assert_relative_eq!(back.arg_periapsis, el.arg_periapsis, max_relative = 1e-9);
        assert_relative_eq!(back.true_anomaly, el.true_anomaly, max_relative = 1e-9);
    }

    #[test]
    fn hyperbolic_round_trip_and_v_infinity() {
        let mu = 1.0;
        let el = Elements { a: -2.0, e: 1.5, inclination: 0.3, raan: 0.2, arg_periapsis: 0.7, true_anomaly: -1.0 };
        let (r, v) = state_from_elements(mu, &el).unwrap();
        let back = elements_from_state(mu, r, v);
        assert_eq!(back.kind, OrbitKind::Hyperbolic);
        assert_relative_eq!(back.elements.a, -2.0, max_relative = 1e-9);
        assert_relative_eq!(back.elements.e, 1.5, max_relative = 1e-9);
        assert_relative_eq!(back.v_infinity.unwrap(), (mu / 2.0f64).sqrt(), max_relative = 1e-9);
        assert_relative_eq!(back.elements.true_anomaly, 2.0 * PI - 1.0, max_relative = 1e-9);
    }

    #[test]
    fn kepler_equation_solvers() {
        for &e in &[0.0, 0.1, 0.5, 0.9, 0.99, 0.999] {
            for k in 0..20 {
                let m = k as f64 * 0.33;
                let ecc = solve_kepler_elliptic(m, e);
                assert!((ecc - e * ecc.sin() - m.rem_euclid(2.0 * PI)).abs() < 1e-12, "e={e} M={m}");
            }
        }
        for &e in &[1.1, 1.5, 3.0, 10.0] {
            for k in -10..10 {
                let m = k as f64 * 1.7;
                let h = solve_kepler_hyperbolic(m, e);
                assert!((e * h.sinh() - h - m).abs() < 1e-9 * m.abs().max(1.0), "e={e} M={m} H={h}");
            }
        }
    }

    #[test]
    fn rutherford_deflection_limits() {
        let mu = 1.0;
        // Impact parameter equal to a = μ/v² gives 90° deflection.
        assert_relative_eq!(hyperbolic_deflection(mu, 1.0, 1.0), PI / 2.0, max_relative = 1e-12);
        assert!(hyperbolic_deflection(mu, 1.0, 1e6) < 1e-5);
        assert_relative_eq!(hyperbolic_periapsis(mu, 1.0, 0.0), 0.0, epsilon = 1e-15);
        let q = hyperbolic_periapsis(mu, 2.0, 3.0);
        assert_relative_eq!(impact_parameter_for_periapsis(mu, 2.0, q), 3.0, max_relative = 1e-12);
    }

    #[test]
    fn scale_lengths() {
        use crate::units::*;
        let hill = hill_radius(AU, 0.0167, M_EARTH, M_SUN);
        assert!((hill / 1.47e9 - 1.0).abs() < 0.02, "Earth Hill radius {hill}");
        let roche = roche_limit_rigid(R_EARTH, 5514.0, 3344.0);
        assert!((roche / 9.5e6 - 1.0).abs() < 0.02, "Moon rigid Roche limit {roche}");
        let p = kepler_period(G_SI * (M_SUN + M_EARTH), AU);
        assert!((p / YEAR - 1.0).abs() < 1e-3);
    }
}
