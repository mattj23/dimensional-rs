use crate::geometry::distances2::dist;
use crate::geometry::line2::intersect_rays;
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

    pub fn spanning_ray(&self, ray: &Ray<f64>) -> Option<Ray<f64>> {
        let results = self.naive_intersections(ray);
        if results.len() == 2 {
            Some(Ray::new(results[0], results[1] - results[0]))

        } else {
            None
        }
    }

    pub fn ray_intersect_with_edge(&self, ray: &Ray<f64>, edge_index: usize) -> Option<f64> {
        let v0 = self.line.points()[self.line.edges()[edge_index].indices.x];
        let v1 = self.line.points()[self.line.edges()[edge_index].indices.y];
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

    pub fn intersect_with_edge(&self, ray: &Ray<f64>, edge_index: usize) -> Option<Point2<f64>> {
        if let Some(t) = self.ray_intersect_with_edge(&ray, edge_index) {
            Some(ray.point_at(t))
        } else {
            None
        }
    }

    pub fn naive_intersections(&self, ray: &Ray<f64>) -> Vec<Point2<f64>> {
        let mut results: Vec<Point2<f64>> = Vec::new();
        for (i, _) in self.line.edges().iter().enumerate() {
            if let Some(point) = self.intersect_with_edge(ray, i) {
                results.push(point);
            }
        }

        results
    }

    // pub fn geometry(&self, ray: &Ray<f64>) -> Vec<Point2<f64>> {
    //     geometry(&self.line, &ray)
    // }
}
