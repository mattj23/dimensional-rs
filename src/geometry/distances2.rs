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

}
