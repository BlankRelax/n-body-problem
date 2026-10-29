//! Barnes–Hut octree (Barnes & Hut, *Nature* 324, 446, 1986).
//!
//! The tree is rebuilt from scratch for every force evaluation (building is
//! O(N log N) and cheap next to the traversal). Nodes are stored in a flat
//! arena; the eight children of a node occupy consecutive slots so a node
//! only stores the index of its first child. Bodies that fall into the same
//! leaf past the maximum depth (coincident or nearly so) are kept in an
//! intrusive linked list and summed directly, so the tree never recurses
//! infinitely.
//!
//! The multipole acceptance criterion is the classic `s / d < θ` with `s`
//! the node width and `d` the distance to the node's centre of mass; only
//! the monopole term is kept. With θ = 0.5 the typical relative force error
//! is ~1e-2, with θ = 0.3 ~1e-3; θ = 0 reproduces direct summation exactly.

use crate::vector::Vec3;

const NONE: u32 = u32::MAX;
const MAX_DEPTH: u32 = 48;

#[derive(Debug, Clone)]
struct Node {
    centre: Vec3,
    half: f64,
    mass: f64,
    /// Mass-weighted position sum during construction, centre of mass after.
    com: Vec3,
    first_child: u32,
    /// Head of the body list for a leaf.
    body: u32,
    count: u32,
}

impl Node {
    fn empty(centre: Vec3, half: f64) -> Node {
        Node { centre, half, mass: 0.0, com: Vec3::ZERO, first_child: NONE, body: NONE, count: 0 }
    }

    #[inline]
    fn is_leaf(&self) -> bool {
        self.first_child == NONE
    }
}

#[derive(Debug, Clone)]
pub struct Octree {
    nodes: Vec<Node>,
    /// `next[b]` links bodies sharing a leaf.
    next: Vec<u32>,
}

#[inline]
fn octant(centre: Vec3, p: Vec3) -> usize {
    (p.x >= centre.x) as usize | ((p.y >= centre.y) as usize) << 1 | ((p.z >= centre.z) as usize) << 2
}

impl Octree {
    pub fn build(pos: &[Vec3], mass: &[f64]) -> Octree {
        let n = pos.len();
        let mut lo = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut hi = -lo;
        for p in pos {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let (centre, half) = if n == 0 {
            (Vec3::ZERO, 1.0)
        } else {
            let c = (lo + hi) * 0.5;
            let extent = (hi - lo).max_abs_component();
            // Slight inflation so that points on the boundary are strictly inside.
            (c, (0.5 * extent * (1.0 + 1e-9)).max(f64::MIN_POSITIVE * 1e10).max(1e-300))
        };
        let mut tree = Octree { nodes: Vec::with_capacity(2 * n + 1), next: vec![NONE; n] };
        tree.nodes.push(Node::empty(centre, half));
        for b in 0..n {
            tree.insert(0, b, 0, pos, mass);
        }
        for node in &mut tree.nodes {
            if node.mass > 0.0 {
                node.com /= node.mass;
            } else if node.count > 0 {
                // Massless bodies: place the node at the first body.
                node.com = pos[node.body as usize];
            }
        }
        tree
    }

    fn subdivide(&mut self, n: usize) -> u32 {
        let first = self.nodes.len() as u32;
        let centre = self.nodes[n].centre;
        let h = self.nodes[n].half * 0.5;
        for o in 0..8 {
            let c = Vec3::new(
                centre.x + if o & 1 != 0 { h } else { -h },
                centre.y + if o & 2 != 0 { h } else { -h },
                centre.z + if o & 4 != 0 { h } else { -h },
            );
            self.nodes.push(Node::empty(c, h));
        }
        self.nodes[n].first_child = first;
        first
    }

    fn insert(&mut self, n: usize, b: usize, depth: u32, pos: &[Vec3], mass: &[f64]) {
        let p = pos[b];
        let m = mass[b];
        {
            let node = &mut self.nodes[n];
            node.mass += m;
            node.com += p * m;
            node.count += 1;
        }
        if self.nodes[n].is_leaf() {
            if self.nodes[n].count == 1 {
                self.nodes[n].body = b as u32;
                self.next[b] = NONE;
                return;
            }
            if depth >= MAX_DEPTH {
                self.next[b] = self.nodes[n].body;
                self.nodes[n].body = b as u32;
                return;
            }
            let first = self.subdivide(n) as usize;
            let centre = self.nodes[n].centre;
            let mut cur = self.nodes[n].body;
            self.nodes[n].body = NONE;
            while cur != NONE {
                let nxt = self.next[cur as usize];
                let o = octant(centre, pos[cur as usize]);
                self.insert(first + o, cur as usize, depth + 1, pos, mass);
                cur = nxt;
            }
            let o = octant(centre, p);
            self.insert(first + o, b, depth + 1, pos, mass);
        } else {
            let centre = self.nodes[n].centre;
            let first = self.nodes[n].first_child as usize;
            self.insert(first + octant(centre, p), b, depth + 1, pos, mass);
        }
    }

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Softened acceleration on body `i` with opening angle `theta`.
    pub fn acceleration_on(&self, i: usize, pos: &[Vec3], mass: &[f64], g: f64, softening: f64, theta: f64, stack: &mut Vec<u32>) -> Vec3 {
        let pi = pos[i];
        let eps2 = softening * softening;
        let mut acc = Vec3::ZERO;
        stack.clear();
        stack.push(0);
        while let Some(n) = stack.pop() {
            let node = &self.nodes[n as usize];
            if node.count == 0 {
                continue;
            }
            if node.is_leaf() {
                let mut cur = node.body;
                while cur != NONE {
                    let j = cur as usize;
                    if j != i {
                        let d = pos[j] - pi;
                        let r2 = d.norm_sq() + eps2;
                        if r2 > 0.0 {
                            acc += d * (g * mass[j] / (r2 * r2.sqrt()));
                        }
                    }
                    cur = self.next[j];
                }
                continue;
            }
            let d = node.com - pi;
            let dist2 = d.norm_sq();
            let s = 2.0 * node.half;
            if s * s < theta * theta * dist2 {
                let r2 = dist2 + eps2;
                acc += d * (g * node.mass / (r2 * r2.sqrt()));
            } else {
                for o in 0..8 {
                    stack.push(node.first_child + o);
                }
            }
        }
        acc
    }
}

/// Barnes–Hut accelerations of every body.
pub fn accelerations(pos: &[Vec3], mass: &[f64], g: f64, softening: f64, theta: f64, out: &mut Vec<Vec3>) {
    let n = pos.len();
    out.clear();
    out.resize(n, Vec3::ZERO);
    if n < 2 {
        return;
    }
    let tree = Octree::build(pos, mass);
    let mut stack = Vec::with_capacity(64);
    for i in 0..n {
        out[i] = tree.acceleration_on(i, pos, mass, g, softening, theta, &mut stack);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::force::direct;
    use rand::{Rng, SeedableRng};

    fn cloud(n: usize, seed: u64) -> (Vec<Vec3>, Vec<f64>) {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let pos = (0..n).map(|_| Vec3::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0))).collect();
        let mass = (0..n).map(|_| rng.gen_range(0.1..1.0)).collect();
        (pos, mass)
    }

    #[test]
    fn theta_zero_matches_direct_summation() {
        let (pos, mass) = cloud(300, 1);
        let (mut a_bh, mut a_d) = (Vec::new(), Vec::new());
        accelerations(&pos, &mass, 1.0, 0.01, 0.0, &mut a_bh);
        direct::accelerations(&pos, &mass, 1.0, 0.01, &mut a_d);
        for (x, y) in a_bh.iter().zip(&a_d) {
            assert!((*x - *y).norm() < 1e-9 * y.norm().max(1.0));
        }
    }

    #[test]
    fn small_theta_is_accurate() {
        let (pos, mass) = cloud(2000, 2);
        let (mut a_bh, mut a_d) = (Vec::new(), Vec::new());
        accelerations(&pos, &mass, 1.0, 0.02, 0.5, &mut a_bh);
        direct::accelerations(&pos, &mass, 1.0, 0.02, &mut a_d);
        let mut worst: f64 = 0.0;
        let mut rms = 0.0;
        for (x, y) in a_bh.iter().zip(&a_d) {
            let rel = (*x - *y).norm() / y.norm();
            worst = worst.max(rel);
            rms += rel * rel;
        }
        rms = (rms / pos.len() as f64).sqrt();
        assert!(rms < 1e-2, "rms relative error {rms}");
        assert!(worst < 0.2, "worst relative error {worst}");
    }

    #[test]
    fn coincident_bodies_terminate() {
        let pos = vec![Vec3::ZERO; 5];
        let mass = vec![1.0; 5];
        let mut a = Vec::new();
        accelerations(&pos, &mass, 1.0, 0.0, 0.5, &mut a);
        assert!(a.iter().all(|v| v.is_finite()));
    }
}
