use ncollide2d::na::{Isometry2, Point2, Vector2};
use ncollide2d::shape::Ball;

pub struct Circle2 {
    pub center: Point2<f64>,
    pub ball: Ball<f64>,
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

    pub fn point_at_angle(&self, angle: f64) -> Point2<f64> {
        let v = Vector2::new(self.ball.radius, 0.0);
        let t = Isometry2::rotation(angle);
        self.center + (t * v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use test_case::test_case;

    #[test_case((0.0, 0.0, 1.0), 0.0, (1.0, 0.0))]
    #[test_case((0.0, 0.0, 1.0), 90.0, (0.0, 1.0))]
    #[test_case((0.0, 0.0, 1.0), 180.0, (-1.0, 0.0))]
    #[test_case((0.0, 0.0, 1.0), 360.0, (1.0, 0.0))]
    #[test_case((1.0, 1.0, 1.0), 0.0, (2.0, 1.0))]
    #[test_case((1.0, 1.0, 1.0), 90.0, (1.0, 2.0))]
    #[test_case((1.0, 1.0, 1.0), 180.0, (0.0, 1.0))]
    #[test_case((1.0, 1.0, 1.0), 360.0, (2.0, 1.0))]
    fn test_circle_point(c: (f64, f64, f64), a: f64, r: (f64, f64)) {
        let circle = Circle2::new(c.0, c.1, c.2);
        let point = circle.point_at_angle(a * std::f64::consts::PI / 180.0);
        assert_relative_eq!(r.0, point.x, epsilon = 1.0e-10);
        assert_relative_eq!(r.1, point.y, epsilon = 1.0e-10);
    }
}
