use ncollide2d::na::{Point2, RealField, Vector2};
use ncollide2d::shape::ConvexPolygon;

/// Return the distance between two 2D points
pub fn dist<N: RealField + Copy>(a: &Point2<N>, b: &Point2<N>) -> N {
    (a - b).norm()
}

pub fn signed_angle<N: RealField + Copy>(v1: &Vector2<N>, v2: &Vector2<N>) -> N {
    (v1.x * v2.y - v1.y * v2.x).atan2(v1.x * v2.x + v1.y * v2.y)
}

/// Find the indices of the farthest pair of points on a convex polygon
pub fn farthest_pair_indices<N: RealField + Copy>(hull: &ConvexPolygon<N>) -> (usize, usize) {
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

    #[test_case((4.1, 3.5, 0.7, 2.4), 0.580388)]
    #[test_case((0.5, -4.6, -2.4, -4.6), -0.589158)]
    #[test_case((1.5, -1.1, 4.0, 4.5), 1.476903)]
    #[test_case((1.8, -0.7, -0.9, -1.8), -1.663553)]
    #[test_case((-3.0, -2.5, -5.0, -3.7), -0.057668)]
    #[test_case((-3.9, -4.6, 3.2, -4.1), 1.365960)]
    #[test_case((-2.0, 3.2, -4.4, -3.7), 1.711390)]
    #[test_case((-5.0, -3.9, 2.5, 0.6), 2.714711)]
    #[test_case((-3.6, 3.0, -1.6, 2.7), -0.341103)]
    #[test_case((0.3, -0.4, -4.0, 0.5), -2.338652)]
    #[test_case((-1.8, 0.8, -4.7, 4.6), -0.356422)]
    #[test_case((0.5, 1.5, -4.2, 2.1), 1.428899)]
    #[test_case((4.6, 1.6, -3.7, -0.7), 2.993835)]
    #[test_case((4.1, -1.6, -0.2, 2.0), 2.042533)]
    #[test_case((2.6, 1.5, -1.9, -1.8), -2.906493)]
    fn test_signed_angles(v: (f64, f64, f64, f64), a: f64) {
        let v1 = Vector2::new(v.0, v.1);
        let v2 = Vector2::new(v.2, v.3);

        assert_relative_eq!(signed_angle(&v1, &v2), a, epsilon = 1e-6);
    }
}
