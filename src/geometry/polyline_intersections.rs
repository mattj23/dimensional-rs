use crate::geometry::distances2::dist;
use crate::geometry::line2::intersect_rays;
use ncollide2d::na::Point2;
use ncollide2d::query::Ray;
use ncollide2d::shape::{ConvexPolygon, Polyline};
use std::error::Error;

pub fn ray_intersect_with_edge(
    line: &Polyline<f64>,
    ray: &Ray<f64>,
    edge_index: usize,
) -> Option<f64> {
    let v0 = line.points()[line.edges()[edge_index].indices.x];
    let v1 = line.points()[line.edges()[edge_index].indices.y];
    let edge_ray = Ray::new(v0, v1 - v0);
    if let Some((t0, t1)) = intersect_rays(&ray, &edge_ray) {
        if t1 >= 0.0 && t1 <= 1.0 {
            Some(t0)
        } else {
            None
        }
    } else {
        None
    }
}

pub fn naive_ray_intersections(line: &Polyline<f64>, ray: &Ray<f64>) -> Vec<f64> {
    let mut results = Vec::new();
    for (i, _) in line.edges().iter().enumerate() {
        if let Some(point) = ray_intersect_with_edge(line, ray, i) {
            results.push(point);
        }
    }

    results
}

/// Finds the projected distance of the farthest point in the polyline from a ray origin in the
/// ray direction
pub fn farthest_point_direction_distance(line: &Polyline<f64>, ray: &Ray<f64>) -> f64 {
    let mut farthest = f64::MIN;
    let n = ray.dir.normalize();
    for v in line.points().iter() {
        farthest = farthest.max(n.dot(&(v - ray.origin)));
    }

    farthest
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ncollide2d::na::Vector2;
    use test_case::test_case;

    #[test_case((3.1, 4.7, 0.9, 3.5), 1.297781)]
    #[test_case((-3.6, -0.6, 2.4, -4.5), 7.517647)]
    #[test_case((1.5, -1.9, -2.6, 0.1), 9.285442)]
    #[test_case((-0.8, 4.5, 3.5, 4.0), 3.076157)]
    #[test_case((3.7, 3.9, 0.4, 3.4), 2.249194)]
    #[test_case((1.4, 4.7, 4.9, -2.0), 6.112484)]
    #[test_case((-0.3, -2.3, 2.3, -0.9), 5.061102)]
    #[test_case((3.6, 2.9, -3.3, -3.3), 12.657211)]
    #[test_case((-1.8, -2.1, -2.7, 0.7), 6.134232)]
    #[test_case((-3.9, -4.5, 1.9, 1.7), 11.441415)]
    fn test_farthest_dist_direction(a: (f64, f64, f64, f64), d: f64) {
        let r = Ray::new(Point2::new(a.0, a.1), Vector2::new(a.2, a.3));
        let result = farthest_point_direction_distance(&sample_polyline(), &r);
        assert_relative_eq!(d, result, epsilon = 1e-5);
    }

    fn sample_polyline() -> Polyline<f64> {
        let points: Vec<Point2<f64>> = vec![
            Point2::new(5.0, 0.0),
            Point2::new(3.5, 0.9),
            Point2::new(4.0, 2.3),
            Point2::new(3.3, 3.3),
            Point2::new(2.7, 4.7),
            Point2::new(1.7, 6.4),
            Point2::new(0.0, 5.9),
            Point2::new(-1.5, 5.7),
            Point2::new(-3.7, 6.4),
            Point2::new(-5.3, 5.3),
            Point2::new(-6.4, 3.7),
            Point2::new(-7.1, 1.9),
            Point2::new(-7.3, 0.0),
            Point2::new(-7.8, -2.1),
            Point2::new(-6.3, -3.7),
            Point2::new(-5.7, -5.7),
            Point2::new(-3.7, -6.3),
            Point2::new(-1.7, -6.2),
            Point2::new(-0.0, -7.2),
            Point2::new(1.5, -5.6),
            Point2::new(2.4, -4.2),
            Point2::new(3.9, -3.9),
            Point2::new(4.9, -2.9),
            Point2::new(4.9, -1.3),
            Point2::new(5.0, 0.0),
        ];
        Polyline::new(points, None)
    }
}
