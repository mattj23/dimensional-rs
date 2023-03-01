use ncollide2d::na::{Point2, Vector2};
use ncollide2d::shape::Ball;

pub struct Circle2 {
    center: Point2<f64>,
    ball: Ball<f64>,
}

impl Circle2 {
    pub fn new(x: f64, y: f64, r: f64) -> Circle2 {
        Circle2 {
            center: Point2::new(x, y),
            ball: Ball::new(r),
        }
    }

    pub fn from_point(center: Point2<f64>, r: f64) -> Circle2 {
        Circle2 {
            center,
            ball: Ball::new(r),
        }
    }
}
