use ncollide2d::na::{Point2, RealField, Vector2};
use ncollide2d::shape::ConvexPolygon;

pub fn dist<N: RealField + Copy>(a: &Point2<N>, b: &Point2<N>) -> N {
    /// Return the distance between two 2D points
    (a - b).norm()
}

pub fn intersection_param<N: RealField + Copy>(
    a0: &Point2<N>,
    ad: &Vector2<N>,
    b0: &Point2<N>,
    bd: &Vector2<N>,
) -> Option<(N, N)> {
    /// Compute the intersection parameters between two parameterized lines. Will return None if
    /// the two directions are parallel to each other
    let det: N = bd.x * ad.y - bd.y * ad.x;
    if det.abs() < N::from_f64(1e-6).unwrap() {
        return Option::None;
    }

    let dx = b0.x - a0.x;
    let dy = b0.y - a0.y;

    Some(((dy * bd.x - dx * bd.y) / det, (dy * ad.x - dx * ad.y) / det))
}

pub fn farthest_pair_indices<N: RealField + Copy>(hull: &ConvexPolygon<N>) -> (usize, usize) {
    /// Find the indices of the farthest pair of points on a convex polygon
    let mut i0: usize = 0;
    let mut i1: usize = 0;
    let mut dist: N = N::from_f64(0.0).unwrap();
    // TODO: Switch to convex hull rotating calipers algorithm
    for i in 0..hull.points().len() {
        for j in 0..hull.points().len() {
            let d: N = (hull.points()[i] - hull.points()[j]).norm();
            if d > dist {
                dist = d;
                i0 = i;
                i1 = j;
            }
        }
    }

    (i0, i1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use test_case::test_case;

    #[test_case((0.0, 0.0, 0.0, 1.0), (1.0, 1.0, 1.0, 0.0), (0.5, 0.5))]
    fn intersection_parameters(av: (f64, f64, f64, f64), bv: (f64, f64, f64, f64), p: (f64, f64)) {
        let a = Point2::new(av.0, av.1);
        let an = Vector2::new(av.2, av.3);
        let b = Point2::new(bv.0, bv.1);
        let bn = Vector2::new(bv.2, bv.3);

        let (ap, bp) = intersection_param(&a, &an, &b, &bn).unwrap();

        assert_relative_eq!(p.0, ap, epsilon = 1.0e-6);
        assert_relative_eq!(p.1, bp, epsilon = 1.0e-6);
    }
}
