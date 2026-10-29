//! A minimal, allocation-free 3-vector in `f64`.
//!
//! Serialises as a 3-element array (`[x, y, z]`) so scenario files stay terse.

use serde::{Deserialize, Serialize};
use std::fmt;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
#[serde(from = "[f64; 3]", into = "[f64; 3]")]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const ZERO: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 0.0 };
    pub const X: Vec3 = Vec3 { x: 1.0, y: 0.0, z: 0.0 };
    pub const Y: Vec3 = Vec3 { x: 0.0, y: 1.0, z: 0.0 };
    pub const Z: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 1.0 };

    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    #[inline]
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    #[inline]
    pub fn cross(self, o: Vec3) -> Vec3 {
        Vec3::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    #[inline]
    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    #[inline]
    pub fn distance(self, o: Vec3) -> f64 {
        (self - o).norm()
    }

    #[inline]
    pub fn distance_sq(self, o: Vec3) -> f64 {
        (self - o).norm_sq()
    }

    /// Unit vector in the same direction; `None` for the zero vector.
    pub fn try_normalized(self) -> Option<Vec3> {
        let n = self.norm();
        if n > 0.0 && n.is_finite() {
            Some(self / n)
        } else {
            None
        }
    }

    /// Unit vector, or the zero vector if the input has zero length.
    pub fn normalized(self) -> Vec3 {
        self.try_normalized().unwrap_or(Vec3::ZERO)
    }

    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    pub fn lerp(self, o: Vec3, t: f64) -> Vec3 {
        self + (o - self) * t
    }

    /// Component of `self` perpendicular to the unit vector `n`.
    pub fn reject_from_unit(self, n: Vec3) -> Vec3 {
        self - n * self.dot(n)
    }

    /// Any unit vector orthogonal to `self` (which need not be normalised).
    /// Uses the "smallest component" trick, which is numerically robust.
    pub fn any_orthonormal(self) -> Vec3 {
        let a = self.x.abs();
        let b = self.y.abs();
        let c = self.z.abs();
        let helper = if a <= b && a <= c {
            Vec3::X
        } else if b <= c {
            Vec3::Y
        } else {
            Vec3::Z
        };
        self.cross(helper).normalized()
    }

    pub fn max_abs_component(self) -> f64 {
        self.x.abs().max(self.y.abs()).max(self.z.abs())
    }

    pub fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

impl From<[f64; 3]> for Vec3 {
    fn from(a: [f64; 3]) -> Self {
        Vec3::new(a[0], a[1], a[2])
    }
}

impl From<Vec3> for [f64; 3] {
    fn from(v: Vec3) -> Self {
        v.to_array()
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    #[inline]
    fn add(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}
impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, o: Vec3) {
        self.x += o.x;
        self.y += o.y;
        self.z += o.z;
    }
}
impl Sub for Vec3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}
impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, o: Vec3) {
        self.x -= o.x;
        self.y -= o.y;
        self.z -= o.z;
    }
}
impl Mul<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, s: f64) -> Vec3 {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }
}
impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline]
    fn mul(self, v: Vec3) -> Vec3 {
        v * self
    }
}
impl MulAssign<f64> for Vec3 {
    #[inline]
    fn mul_assign(&mut self, s: f64) {
        self.x *= s;
        self.y *= s;
        self.z *= s;
    }
}
impl Div<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn div(self, s: f64) -> Vec3 {
        Vec3::new(self.x / s, self.y / s, self.z / s)
    }
}
impl DivAssign<f64> for Vec3 {
    #[inline]
    fn div_assign(&mut self, s: f64) {
        self.x /= s;
        self.y /= s;
        self.z /= s;
    }
}
impl Neg for Vec3 {
    type Output = Vec3;
    #[inline]
    fn neg(self) -> Vec3 {
        Vec3::new(-self.x, -self.y, -self.z)
    }
}
impl Index<usize> for Vec3 {
    type Output = f64;
    fn index(&self, i: usize) -> &f64 {
        match i {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Vec3 index out of range: {i}"),
        }
    }
}
impl IndexMut<usize> for Vec3 {
    fn index_mut(&mut self, i: usize) -> &mut f64 {
        match i {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Vec3 index out of range: {i}"),
        }
    }
}
impl Sum for Vec3 {
    fn sum<I: Iterator<Item = Vec3>>(iter: I) -> Vec3 {
        iter.fold(Vec3::ZERO, |a, b| a + b)
    }
}
impl<'a> Sum<&'a Vec3> for Vec3 {
    fn sum<I: Iterator<Item = &'a Vec3>>(iter: I) -> Vec3 {
        iter.fold(Vec3::ZERO, |a, b| a + *b)
    }
}
impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(p) = f.precision() {
            write!(f, "({:.*e}, {:.*e}, {:.*e})", p, self.x, p, self.y, p, self.z)
        } else {
            write!(f, "({:.4e}, {:.4e}, {:.4e})", self.x, self.y, self.z)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_product_is_right_handed() {
        assert_eq!(Vec3::X.cross(Vec3::Y), Vec3::Z);
        assert_eq!(Vec3::Y.cross(Vec3::Z), Vec3::X);
        assert_eq!(Vec3::Z.cross(Vec3::X), Vec3::Y);
    }

    #[test]
    fn orthonormal_helper_is_orthogonal_and_unit() {
        for v in [Vec3::new(1.0, 2.0, 3.0), Vec3::X, Vec3::new(0.0, 0.0, -5.0), Vec3::new(1e-9, 1.0, 1e-9)] {
            let o = v.any_orthonormal();
            assert!(o.dot(v).abs() < 1e-12 * v.norm().max(1.0));
            assert!((o.norm() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn serde_round_trip_as_array() {
        let v = Vec3::new(1.5, -2.0, 3.25);
        let s = serde_json::to_string(&v).unwrap();
        assert_eq!(s, "[1.5,-2.0,3.25]");
        let back: Vec3 = serde_json::from_str(&s).unwrap();
        assert_eq!(v, back);
    }
}
