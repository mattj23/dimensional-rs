use crate::intersections::distances2::dist;
use ncollide2d::na::Point2;
use ncollide2d::procedural::Polyline;
use ncollide2d::query::Ray;
use ncollide2d::shape::ConvexPolygon;
use std::error::Error;

pub struct ClosedPolyline {
    pub line: Polyline<f64>,
    pub hull: ConvexPolygon<f64>,
}

impl ClosedPolyline {
    pub fn new(
        points: &Vec<Point2<f64>>,
        tol: Option<f64>,
    ) -> Result<ClosedPolyline, Box<dyn Error>> {
        /// Create a new closed polyline. The provided points will be copied and adjacent duplicates
        /// removed.  If the provided points are not enclosed, the first point will be copied to the
        /// end.
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

    // pub fn intersect_with_edge(&self, ray: &Ray<f64>, edge_index: usize) -> Option<Point2<f64>> {
    //     intersect_with_edge(&self.line, &ray, edge_index)
    // }

    // pub fn naive_intersections(&self, ray: &Ray<f64>) -> Vec<Point2<f64>> {
    //     let mut results: Vec<Point2<f64>> = Vec::new();
    //     for (i, _) in self.line.edges().iter().enumerate() {
    //         if let Some(point) = self.intersect_with_edge(ray, i) {
    //             results.push(point);
    //         }
    //     }
    //
    //     results
    // }

    // pub fn intersections(&self, ray: &Ray<f64>) -> Vec<Point2<f64>> {
    //     intersections(&self.line, &ray)
    // }
}
