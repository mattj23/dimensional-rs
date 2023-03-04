use crate::algorithms::preceding_index_search;
use crate::errors::InvalidGeometry;
use crate::geometry::distances2::{dist, signed_angle};
use ncollide2d::math::Isometry;
use ncollide2d::na::{Point2, Unit, Vector2};
use ncollide2d::shape::Polyline;
use std::error::Error;
use std::ops::Index;

type UnitVec2 = Unit<Vector2<f64>>;

fn sym_unit_vec(a: &UnitVec2, b: &UnitVec2) -> UnitVec2 {
    let t = signed_angle(a, b) * 0.5;
    Isometry::rotation(t) * a
}

/// A Curve2 is a 2 dimensional polygonal chain in which its points are connected. It optionally
/// may include normals. This struct and its methods allow for convenient handling of distance
/// searches, transformations, resampling, and splitting.
pub struct Curve2 {
    pub line: Polyline<f64>,
    normals: Vec<UnitVec2>,
    lengths: Vec<f64>,
    is_closed: bool,
    tol: f64,
}

impl Curve2 {
    fn first_n(&self) -> UnitVec2 {
        self.line.edges().first().unwrap().normal.unwrap()
    }

    fn last_n(&self) -> UnitVec2 {
        self.line.edges().last().unwrap().normal.unwrap()
    }

    fn last_vi(&self) -> usize {
        self.line.points().len() - 1
    }

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
        let d = l.clamp(0.0, self.length());

        // We know that we have at least two vertices, and that the index that will be returned is
        // going to be between 0 and self.line.points.len() - 1. We can check the end cases first.
        let index = preceding_index_search(&self.lengths, l);

        if index == self.line.points().len() - 1 {
            (index - 1, 0.0)
        } else {
            let f = (d - self.lengths[index]) / (self.lengths[index + 1] - self.lengths[index]);
            (index, 1.0 - f)
        }
    }

    pub fn point_at(&self, l: f64) -> Point2<f64> {
        let (i, f) = self.at_length(l);
        let p = self.line.points()[i];
        let v = self.line.points()[i + 1] - p;
        p + (1.0 - f) * v
    }

    fn normal_at_vertex(&self, i: usize) -> UnitVec2 {
        if self.is_closed && (i == 0 || i == self.last_vi()) {
            sym_unit_vec(&self.first_n(), &self.last_n())
        } else if i == 0 {
            self.first_n()
        } else if i == self.last_vi() {
            self.last_n()
        } else {
            let n0 = self.line.edges()[i - 1].normal.unwrap();
            let n1 = self.line.edges()[i].normal.unwrap();

            sym_unit_vec(&n0, &n1)
        }
    }

    pub fn normal_at(&self, l: f64) -> UnitVec2 {
        let (i, f) = self.at_length(l);

        if f <= self.tol {
            self.normal_at_vertex(i)
        } else if (f - 1.0).abs() <= self.tol {
            self.normal_at_vertex(self.last_vi())
        } else {
            self.line.edges()[i].normal.unwrap()
        }
    }

    pub fn point_at_fraction(&self, f: f64) -> Point2<f64> {
        self.point_at(f / self.length())
    }

    pub fn from_points(
        points: &[Point2<f64>],
        tol: f64,
        force_closed: bool,
    ) -> Result<Self, Box<dyn Error>> {
        let mut pts = points.to_vec();
        pts.dedup_by(|a, b| dist(a, b) <= tol);

        if pts.len() < 2 {
            return Err(Box::try_from(InvalidGeometry::NotEnoughPoints).unwrap());
        }

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

        Ok(Curve2 {
            line,
            normals,
            lengths,
            is_closed,
            tol,
        })
    }

    pub fn edge_count(&self) -> usize {
        self.line.edges().len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ncollide2d::na::Point2;
    use test_case::test_case;

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
        let curve = Curve2::from_points(&points, 1e-6, false).unwrap();

        assert!(!curve.is_closed());
        assert_relative_eq!(3.0, curve.length(), epsilon = 1e-10);
    }

    #[test]
    fn test_create_force_closed() {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();

        assert!(curve.is_closed());
        assert_relative_eq!(4.0, curve.length(), epsilon = 1e-10);
    }

    #[test]
    fn test_create_naturally_closed() {
        let points = sample_points(&sample2());
        let curve = Curve2::from_points(&points, 1e-6, false).unwrap();

        assert!(curve.is_closed());
        assert_relative_eq!(4.0, curve.length(), epsilon = 1e-10);
    }

    #[test_case(0.5, 0, 0.5)]
    #[test_case(-0.5, 0, 1.0)]
    #[test_case(0.0, 0, 1.0)]
    #[test_case(5.0, 3, 0.0)]
    #[test_case(2.0, 2, 1.0)]
    #[test_case(2.25, 2, 0.75)]
    fn test_lengths(l: f64, ei: usize, ef: f64) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();

        let (i, f) = curve.at_length(l);
        assert_eq!(ei, i);
        assert_relative_eq!(ef, f, epsilon = 1e-8);
    }

    #[test_case(0.5, (0.5, 0.0))]
    #[test_case(-0.5, (0.0, 0.0))]
    #[test_case(0.0, (0.0, 0.0))]
    #[test_case(5.0, (0.0, 0.0))]
    #[test_case(2.0, (1.0, 1.0))]
    #[test_case(2.25, (0.75, 1.0))]
    fn test_points_at_length(l: f64, e: (f64, f64)) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();
        let result = curve.point_at(l);

        assert_relative_eq!(e.0, result.x, epsilon = 1e-8);
        assert_relative_eq!(e.1, result.y, epsilon = 1e-8);
    }
}
