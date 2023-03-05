use crate::geometry::distances2::signed_angle;
use ncollide2d::math::Isometry;
use ncollide2d::na::{Unit, Vector2};
use std::cmp::Ordering;

// TODO: Move these to a common module?
pub type UnitVec2 = Unit<Vector2<f64>>;

pub fn sym_unit_vec(a: &UnitVec2, b: &UnitVec2) -> UnitVec2 {
    let t = signed_angle(a, b) * 0.5;
    Isometry::rotation(t) * a
}

#[derive(Copy, Clone)]
pub struct IndAndFrac {
    pub i: usize,
    pub f: f64,
}

impl IndAndFrac {
    pub fn new(i: usize, f: f64) -> Self {
        Self {
            i,
            f: f.clamp(0.0, 1.0),
        }
    }

    pub fn one(i: usize) -> Self {
        Self { i, f: 1.0 }
    }

    pub fn zero(i: usize) -> Self {
        Self { i, f: 0.0 }
    }

    pub fn next_zero(&self) -> Self {
        Self::zero(self.i + 1)
    }
}

impl Eq for IndAndFrac {}

impl PartialEq<Self> for IndAndFrac {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.f == other.f
    }
}

impl PartialOrd<Self> for IndAndFrac {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if other.i != self.i {
            Some(self.i.cmp(&other.i))
        } else {
            self.f.partial_cmp(&other.f)
        }
    }
}

impl Ord for IndAndFrac {
    fn cmp(&self, other: &Self) -> Ordering {
        if other.i != self.i {
            self.i.cmp(&other.i)
        } else {
            self.f.partial_cmp(&other.f).unwrap()
        }
    }
}
