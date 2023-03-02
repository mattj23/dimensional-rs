use ncollide2d::na::{Isometry2, Point2, RealField, Vector2};
use ncollide2d::query::Ray;
use ncollide2d::shape::Ball;

/// Compute the intersection parameters between two parameterized lines. Will return None if
/// the two directions are parallel to each other
pub fn intersection_param<N: RealField + Copy>(
    a0: &Point2<N>,
    ad: &Vector2<N>,
    b0: &Point2<N>,
    bd: &Vector2<N>,
) -> Option<(N, N)> {
    let det: N = bd.x * ad.y - bd.y * ad.x;
    if det.abs() < N::from_f64(1e-6).unwrap() {
        return None;
    }

    let dx = b0.x - a0.x;
    let dy = b0.y - a0.y;

    Some(((dy * bd.x - dx * bd.y) / det, (dy * ad.x - dx * ad.y) / det))
}

pub fn intersect_rays(r0: &Ray<f64>, r1: &Ray<f64>) -> Option<(f64, f64)> {
    intersection_param(&r0.origin, &r0.dir, &r1.origin, &r1.dir)
}

pub trait Line2 {
    fn origin(&self) -> Point2<f64>;
    fn dir(&self) -> Vector2<f64>;
    fn at(&self, t: f64) -> Point2<f64>;

    fn projected_parameter(&self, p: &Point2<f64>) -> f64 {
        let n = p - self.origin();
        self.dir().dot(&n)
    }

    fn projected_point(&self, p: &Point2<f64>) -> Point2<f64> {
        self.at(self.projected_parameter(p))
    }

    /// Returns the direction of the vector turned in its orthogonal direction by rotating it
    /// 90 degrees
    fn orthogonal(&self) -> Vector2<f64> {
        Isometry2::rotation(std::f64::consts::PI / 2.0) * self.dir()
    }

    /// Returns a ray that has rotated this entity by 90 degrees
    fn turned(&self) -> Ray<f64> {
        Ray::new(self.origin(), self.orthogonal())
    }

    // fn intersect_line(&self, other: &impl Line2) -> Option<(f64, f64)> {
    //     intersection_param(&self.origin(), &self.dir(), &other.origin(), &other.dir())
    // }
}

impl Line2 for Ray<f64> {
    fn origin(&self) -> Point2<f64> {
        self.origin
    }

    fn dir(&self) -> Vector2<f64> {
        self.dir
    }

    fn at(&self, t: f64) -> Point2<f64> {
        self.point_at(t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use test_case::test_case;

    /// These tests check that the intersection parameter calculation between two parameterized
    /// lines works as expected. The test cases were generated by starting with a random
    /// intersection point and selecting random orthogonal vector and parameter values that rounded
    /// to a single decimal place.
    #[test_case((11.0, 0.7, -4.2, -2.7), (-0.1, -4.7, 1.8, 0.0), (2.0, 1.5))]
    #[test_case((6.2, 3.0, 1.3, -1.9), (-1.3, 8.3, 3.1, -1.7), (-1.0, 2.0))]
    #[test_case((3.6, -3.3, -4.1, 3.2), (9.9, 0.9, 3.0, 2.0), (0.0, -2.1))]
    #[test_case((7.4, -10.1, -2.4, 3.2), (-7.0, -7.7, 4.2, 2.8), (2.5, 2.0))]
    #[test_case((-1.6, 4.3, -0.0, -2.3), (-9.8, 7.1, -4.1, 3.7), (2.0, -2.0))]
    #[test_case((1.7, -0.5, -0.5, -3.3), (3.9, -2.3, 0.5, -1.5), (-1.0, -3.4))]
    #[test_case((4.7, 0.4, -0.6, 0.9), (4.7, 1.3, -3.0, -0.0), (1.0, 0.2))]
    #[test_case((-4.8, -2.0, -1.1, -1.3), (-0.8, -21.0, -0.4, 4.8), (-2.0, 4.5))]
    #[test_case((9.1, 5.7, 4.9, 1.6), (-9.1, -15.5, 2.1, 4.5), (-2.0, 4.0))]
    #[test_case((2.8, 15.7, 0.7, 3.0), (-3.9, -7.1, 1.3, 3.6), (-4.0, 3.0))]
    #[test_case((-5.0, 6.4, -2.6, 0.6), (5.3, 2.3, -2.4, 4.0), (-3.5, 0.5))]
    #[test_case((0.4, -1.9, 2.5, 1.5), (7.6, -19.0, -3.2, 4.2), (-1.6, 3.5))]
    #[test_case((10.6, 5.9, 2.0, 0.7), (-2.0, 7.1, 2.6, -4.7), (-5.0, 1.0))]
    #[test_case((7.3, -0.1, -3.5, 0.9), (14.4, 14.6, -4.4, -3.0), (3.0, 4.0))]
    #[test_case((-5.3, 11.5, -0.8, 4.6), (-15.7, 10.0, -3.1, 2.5), (-2.5, -4.0))]
    #[test_case((4.8, 1.9, 1.0, 1.0), (4.8, 1.9, -1.1, 1.0), (0.0, 0.0))]
    fn inter_param_success(av: (f64, f64, f64, f64), bv: (f64, f64, f64, f64), p: (f64, f64)) {
        let a = Point2::new(av.0, av.1);
        let an = Vector2::new(av.2, av.3);
        let b = Point2::new(bv.0, bv.1);
        let bn = Vector2::new(bv.2, bv.3);

        let (ap, bp) = intersection_param(&a, &an, &b, &bn).unwrap();

        assert_relative_eq!(p.0, ap, epsilon = 1.0e-6);
        assert_relative_eq!(p.1, bp, epsilon = 1.0e-6);
    }

    /// These tests check that the intersection parameter calculation between two parameterized
    /// lines which are parallel returns a None
    #[test_case((-5.0, 2.8, 2.2, 1.8), (-4.2, -0.2, 6.6, 5.4))]
    #[test_case((3.3, 2.5, 4.0, 1.0), (3.2, 0.7, -20.0, -5.0))]
    #[test_case((4.2, -2.3, -0.6, 1.4), (-1.0, 0.5, -2.4, 5.6))]
    #[test_case((-1.1, 2.0, 5.0, 4.0), (4.9, -2.8, 19.5, 15.6))]
    #[test_case((2.4, -3.0, -1.8, -2.6), (0.1, 0.7, 7.2, 10.4))]
    #[test_case((1.2, 2.1, 4.3, -1.0), (-1.4, 3.9, 8.6, -2.0))]
    #[test_case((-4.8, -2.0, -0.1, -0.5), (3.0, -3.6, 0.4, 2.0))]
    #[test_case((-4.4, -0.4, 3.1, 1.1), (3.1, 4.9, -3.1, -1.1))]
    #[test_case((-1.0, -0.1, 0.1, 1.0), (1.2, -2.8, 0.3, 3.0))]
    #[test_case((4.7, -3.7, 4.0, -2.5), (2.5, -0.4, -11.2, 7.0))]
    #[test_case((-1.2, 0.4, 3.9, 0.9), (2.7, -4.6, -11.7, -2.7))]
    #[test_case((-4.8, 4.1, 4.4, -3.0), (4.3, 3.2, -6.6, 4.5))]
    #[test_case((-3.7, 4.7, -1.4, 2.6), (1.0, -4.5, -2.1, 3.9))]
    #[test_case((-0.1, 2.6, 1.4, -0.3), (-0.6, -1.5, -4.2, 0.9))]
    #[test_case((1.8, -2.0, 4.5, -2.0), (0.5, -4.7, 18.0, -8.0))]
    fn inter_parallel_fail(av: (f64, f64, f64, f64), bv: (f64, f64, f64, f64)) {
        let a = Point2::new(av.0, av.1);
        let an = Vector2::new(av.2, av.3);
        let b = Point2::new(bv.0, bv.1);
        let bn = Vector2::new(bv.2, bv.3);

        let result = intersection_param(&a, &an, &b, &bn);

        assert_eq!(None, result);
    }
}
