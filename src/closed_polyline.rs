use crate::geometry::distances2::dist;
use crate::geometry::line2::intersect_rays;
use crate::geometry::polyline_intersections::naive_ray_intersections;
use ncollide2d::na::Point2;
use ncollide2d::query::Ray;
use ncollide2d::shape::{ConvexPolygon, Polyline};
use std::error::Error;

pub struct ClosedPolyline {
    pub line: Polyline<f64>,
    pub hull: ConvexPolygon<f64>,
}

impl ClosedPolyline {
    /// Create a new closed polyline. The provided points will be copied and adjacent duplicates
    /// removed.  If the provided points are not enclosed, the first point will be copied to the
    /// end.
    pub fn new(
        points: &Vec<Point2<f64>>,
        tol: Option<f64>,
    ) -> Result<ClosedPolyline, Box<dyn Error>> {
        // Construct the vertices we're looking for
        let mut vertices: Vec<Point2<f64>> = vec![points[0]];
        let tol_value = tol.unwrap_or(1e-6);

        for i in 0..points.len() - 1 {
            if dist(&points[i], &points[i + 1]) >= tol_value {
                vertices.push(points[i]);
            }
        }

        if dist(&vertices[0], vertices.last().unwrap()) > tol_value {
            vertices.push(vertices[0]);
        }

        let hull: ConvexPolygon<f64> = ConvexPolygon::try_from_points(&vertices).unwrap();
        let line: Polyline<f64> = Polyline::new(vertices, Option::None);

        Ok(ClosedPolyline { line, hull })
    }

    /// Attempts to create a "spanning ray" through the polyline along the parameterized line
    /// represented by the ray argument. A "spanning ray" is a ray that starts on the surface of
    /// the polyline and passes through it ending at the other side, such that t=0 is an
    /// intersection with the polyline on one side, t=1.0 is an intersection on the other side, and
    /// there are no additional intersections between them. The spanning ray will have the same
    /// direction as the original intersection ray.
    pub fn spanning_ray(&self, ray: &Ray<f64>) -> Option<Ray<f64>> {
        let mut results = naive_ray_intersections(&self.line, ray);
        results.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if results.len() == 2 {
            let p0 = ray.point_at(results[0]);
            let p1 = ray.point_at(results[1]);
            Some(Ray::new(p0, p1 - p0))
        } else {
            None
        }
    }
}
