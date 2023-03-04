use crate::geometry::distances2::dist;
use ncollide2d::na::{Point2, Unit, Vector2};
use ncollide2d::shape::Polyline;

type UnitVec2 = Unit<Vector2<f64>>;

/// A Curve2 is a 2 dimensional polygonal chain in which its points are connected. It optionally
/// may include normals. This struct and its methods allow for convenient handling of distance
/// searches, transformations, resampling, and splitting.
pub struct Curve2 {
    pub line: Polyline<f64>,
    normals: Vec<UnitVec2>,
    lengths: Vec<f64>,
    is_closed: bool,
    tol: f64
}

impl Curve2 {
    pub fn is_closed(&self) -> bool {
        self.is_closed
    }

    pub fn length(&self) -> f64 {
        *self.lengths.last().unwrap_or(&0.0)
    }

    /// Finds the preceding index of the vertex at the given length along the curve, and returns
    /// an f64 which represents the weighting of the preceding point needed to reconstruct the
    /// properties of the point at the specified length.
    ///
    /// For example, imagine a curve of length 1 which has 3 vertices. At l<=0.0, the function will
    /// return (0, 1.0). At l>=1.0, it will return (2, 0.0).  The position of the point at l can
    /// be found by self.line.points[i] * f + self.line.points[i+1] * (1.0 - f), where i is the
    /// index and f is the f64
    fn at_length(&self, l: f64) -> (usize, f64) {
        todo!()
    }

    pub fn from_points(points: &[Point2<f64>], tol: f64, force_closed: bool) -> Self {
        let mut pts = points.to_vec();
        pts.dedup_by(|a, b| dist(a, b) <= tol);

        // Check if the curve is supposed to be closed
        if let (true, Some(start), Some(end)) = (force_closed, pts.first(), pts.last()) {
            if dist(start, end) > tol {
                pts.push(*start);
            }
        }

        let is_closed = pts.len() >= 2 && dist(&pts[0], pts.last().unwrap()) <= tol;

        // Because we will not actually pass indices into the polyline creation method we can
        // trust that the edges will match the vertex indices.  There will be one less edge than
        // there is vertices, and each edge i will join vertex i with vertex i+1
        let line = Polyline::new(pts, None);

        // Compute the normals
        let normals: Vec<UnitVec2> = line.edges().iter().map(|e| e.normal.unwrap()).collect();

        let mut lengths: Vec<f64> = vec![0.0];
        for e in line.edges().iter() {
            let d = dist(&line.points()[e.indices[0]], &line.points()[e.indices[1]]);
            lengths.push(d + lengths.last().unwrap_or(&0.0));
        }

        Curve2 {
            line,
            normals,
            lengths,
            is_closed,
            tol,
        }
    }

    pub fn edge_count(&self) -> usize {
        self.line.edges().len()
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use super::*;
    use ncollide2d::na::Point2;

    fn sample1() -> Vec<(f64, f64)> {
        vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    }

    fn sample2() -> Vec<(f64, f64)> {
        vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
    }


    fn sample_points(p: &[(f64, f64)]) -> Vec<Point2<f64>> {
        p.iter().map(|(a, b)| Point2::new(*a, *b)).collect()
    }

    #[test]
    fn test_create_open() {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, false);

        assert!(!curve.is_closed());
        assert_relative_eq!(3.0, curve.length(), epsilon = 1e-10);
    }

    #[test]
    fn test_create_force_closed() {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true);

        assert!(curve.is_closed());
        assert_relative_eq!(4.0, curve.length(), epsilon = 1e-10);
    }

    #[test]
    fn test_create_naturally_closed() {
        let points = sample_points(&sample2());
        let curve = Curve2::from_points(&points, 1e-6, false);

        assert!(curve.is_closed());
        assert_relative_eq!(4.0, curve.length(), epsilon = 1e-10);
    }
}
