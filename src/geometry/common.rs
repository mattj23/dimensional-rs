use crate::geometry::distances2::signed_angle;
use ncollide2d::math::Isometry;
use ncollide2d::na::{Unit, Vector2};

// TODO: Move these to a common module?
pub type UnitVec2 = Unit<Vector2<f64>>;

pub fn sym_unit_vec(a: &UnitVec2, b: &UnitVec2) -> UnitVec2 {
    let t = signed_angle(a, b) * 0.5;
    Isometry::rotation(t) * a
}

pub struct IndexAndFraction {
    pub i: usize,
    pub f: f64,
}

impl IndexAndFraction {
    pub fn new(i: usize, f: f64) -> Self {
        Self { i, f }
    }
}
