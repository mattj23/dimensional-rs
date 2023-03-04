use crate::geometry::distances2::{dist, mid_point, signed_angle};
use crate::geometry::line2::{intersect_rays, Line2};
use crate::serialize::Ray2f64;
use ncollide2d::na::{Isometry2, Point2, Vector2};
use ncollide2d::query::Ray;
use ncollide2d::shape::Polyline;
use serde::Serialize;

pub fn cleaned_polyline(points: &[Point2<f64>], tol: f64) -> Polyline<f64> {
    let mut vertices = points.to_vec();
    vertices.dedup_by(|a, b| dist(a, b) <= tol);
    Polyline::new(vertices, Option::None)
}

#[derive(Clone, Serialize)]
pub struct SpanningRay {
    #[serde(with = "Ray2f64")]
    ray: Ray<f64>,
}

impl SpanningRay {
    pub fn new(p0: Point2<f64>, p1: Point2<f64>) -> SpanningRay {
        SpanningRay {
            ray: Ray::new(p0, p1 - p0),
        }
    }

    pub fn symmetry(&self, other: &SpanningRay) -> Ray<f64> {
        let angle = signed_angle(&self.ray.dir, &other.ray.dir) * 0.5;
        Ray::new(
            mid_point(&self.ray.origin, &other.ray.origin),
            Isometry2::rotation(angle) * self.ray.dir,
        )
    }

    pub fn reversed(&self) -> SpanningRay {
        SpanningRay::new(self.ray.point_at(1.0), self.ray.origin)
    }
}

impl Line2 for SpanningRay {
    fn origin(&self) -> Point2<f64> {
        self.ray.origin
    }

    fn dir(&self) -> Vector2<f64> {
        self.ray.dir
    }

    fn at(&self, t: f64) -> Point2<f64> {
        self.ray.point_at(t)
    }
}

pub fn ray_intersect_with_edge(
    line: &Polyline<f64>,
    ray: &Ray<f64>,
    edge_index: usize,
) -> Option<f64> {
    let v0 = line.points()[line.edges()[edge_index].indices.x];
    let v1 = line.points()[line.edges()[edge_index].indices.y];
    let edge_ray = Ray::new(v0, v1 - v0);
    if let Some((t0, t1)) = intersect_rays(ray, &edge_ray) {
        if (0.0..=1.0).contains(&t1) {
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

pub fn max_intersection(line: &Polyline<f64>, ray: &Ray<f64>) -> Option<f64> {
    let ts = naive_ray_intersections(line, ray);
    ts.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).cloned()
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

/// Attempts to create a "spanning ray" through the polyline along the parameterized line
/// represented by the ray argument. A "spanning ray" is a ray that starts on the surface of
/// the polyline and passes through it ending at the other side, such that t=0 is an
/// intersection with the polyline on one side, t=1.0 is an intersection on the other side, and
/// there are no additional intersections between them. The spanning ray will have the same
/// direction as the original intersection ray.
pub fn spanning_ray(line: &Polyline<f64>, ray: &Ray<f64>) -> Option<SpanningRay> {
    let mut results = naive_ray_intersections(line, ray);
    results.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if results.len() == 2 {
        Some(SpanningRay::new(
            ray.point_at(results[0]),
            ray.point_at(results[1]),
        ))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ncollide2d::na::{Point2, Vector2};
    use test_case::test_case;

    #[test]
    fn test_symmetry_ray() {
        let s0 = SpanningRay::new(Point2::new(1.0, 0.0), Point2::new(2.0, 0.0));
        let s1 = SpanningRay::new(Point2::new(0.0, 1.0), Point2::new(0.0, 2.0));
        let r = s0.symmetry(&s1);
        let en = Vector2::new(1.0, 1.0).normalize();

        assert_relative_eq!(0.5, r.origin.x, epsilon = 1e-8);
        assert_relative_eq!(0.5, r.origin.y, epsilon = 1e-8);

        assert_relative_eq!(en.x, r.dir.x, epsilon = 1e-8);
        assert_relative_eq!(en.y, r.dir.y, epsilon = 1e-8);
    }

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
