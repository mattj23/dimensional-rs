use crate::algorithms::preceding_index_search;
use crate::errors::InvalidGeometry;
use crate::geometry::common::{sym_unit_vec, IndAndFrac, UnitVec2};
use crate::geometry::distances2::dist;
use crate::geometry::line2::Line2;
use crate::geometry::partitioning::{BreadthFirst, PointVisitor, SearchType};
use crate::geometry::polyline::{max_intersection, spanning_ray, SpanningRay};
use itertools::Itertools;
use ncollide2d::math::Isometry;
use ncollide2d::na::Point2;
use ncollide2d::query::Ray;
use ncollide2d::shape::{ConvexPolygon, Polyline};
use std::error::Error;

/// A Curve2 is a 2 dimensional polygonal chain in which its points are connected. It optionally
/// may include normals. This struct and its methods allow for convenient handling of distance
/// searches, transformations, resampling, and splitting.
pub struct Curve2 {
    pub line: Polyline<f64>,
    lengths: Vec<f64>,
    is_closed: bool,
    tol: f64,
}

impl Curve2 {
    pub fn lengths(&self) -> &Vec<f64> {
        &self.lengths
    }

    pub fn tol(&self) -> f64 {
        self.tol
    }
    pub fn first_n(&self) -> UnitVec2 {
        self.line.edges().first().unwrap().normal.unwrap()
    }

    pub fn last_n(&self) -> UnitVec2 {
        self.line.edges().last().unwrap().normal.unwrap()
    }

    pub fn first_v(&self) -> Point2<f64> {
        self.line.points()[0]
    }

    pub fn last_v(&self) -> Point2<f64> {
        *self.line.points().last().unwrap()
    }

    pub fn first_d(&self) -> UnitVec2 {
        dir_from_normal(&self.first_n())
    }

    pub fn last_d(&self) -> UnitVec2 {
        dir_from_normal(&self.last_n())
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

    fn at_length(&self, l: f64) -> IndAndFrac {
        let d = l.clamp(0.0, self.length());

        // We know that we have at least two vertices, and that the index that will be returned is
        // going to be between 0 and self.line.points.len() - 1. We can check the end cases first.
        let index = preceding_index_search(&self.lengths, l);

        if index == self.line.points().len() - 1 {
            IndAndFrac::one(index - 1)
        } else {
            let f = (d - self.lengths[index]) / (self.lengths[index + 1] - self.lengths[index]);
            IndAndFrac::new(index, f)
        }
    }

    fn point_from_index_fraction(&self, iaf: &IndAndFrac) -> Point2<f64> {
        let p = self.line.points()[iaf.i];
        let v = self.line.points()[iaf.i + 1] - p;
        p + iaf.f * v
    }

    pub fn point_at(&self, l: f64) -> Point2<f64> {
        self.point_from_index_fraction(&self.at_length(l))
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

    /// Finds the closest point on the curve to the query point, returning the index of the
    /// edge and the location of the found point.
    fn closest_point_and_edge(&self, query: &Point2<f64>) -> (usize, Point2<f64>) {
        let mut collected = Vec::new();
        let mut visitor = PointVisitor::new(query, &mut collected, SearchType::Closest);
        self.line.bvt().bf_visit(&mut visitor);

        let mut result: Option<(f64, usize, Point2<f64>)> = None;
        for t in collected.iter() {
            let e = &self.line.edges()[t.value];
            let p0 = self.line.points()[e.indices.coords.x];
            let p1 = self.line.points()[e.indices.coords.y];
            let r = Ray::new(p0, p1 - p0);
            let r_t = r.projected_parameter(query);
            let cp = r.point_at(r_t.clamp(0.0, 1.0));
            let d = dist(&cp, query);
            if let Some(value) = result {
                if d < value.0 {
                    result = Some((d, t.value, cp));
                }
            } else {
                result = Some((d, t.value, cp));
            }
        }

        let value = result.unwrap();
        (value.1, value.2)
    }

    pub fn closest_point_to(&self, query: &Point2<f64>) -> Point2<f64> {
        self.closest_point_and_edge(query).1
    }

    /// Projects the point onto the curve and returns the distance along the curve that the
    /// projected point is
    pub fn distance_along(&self, point: &Point2<f64>) -> f64 {
        let (index, p) = self.closest_point_and_edge(point);
        self.lengths[index] + dist(&p, &self.line.points()[index])
    }

    /// Returns a curve portion between the section at length l0 and l1. If the curve is not closed,
    /// the case where l1 < l0 will return None. If the curve is closed, the portion of the curve
    /// which is returned will depend on whether l0 is larger or smaller than l1.
    ///
    /// The new curve will begin at the point corresponding with l0.
    pub fn portion_between_lengths(&self, l0: f64, l1: f64) -> Option<Curve2> {
        // If either the distance between l1 and l0 are less than the curve tolerance or the orders
        // are inverted when the curve isn't closed, we have a poorly conditioned request and we
        // can return None

        let start = self.at_length(l0);
        let end = self.at_length(l1);
        let mut wrap = end < start;

        if (l1 - l0).abs() < self.tol || (!self.is_closed && wrap) {
            None
        } else {
            let mut points = Vec::new();
            let mut working = start;

            loop {
                points.push(self.point_from_index_fraction(&working));
                working = working.next_zero();

                // Terminal condition
                if working >= end && !wrap {
                    break;
                }

                if working.i >= self.last_vi() && wrap {
                    working = IndAndFrac::zero(0);
                    wrap = false;
                }
            }
            points.push(self.point_from_index_fraction(&end));

            if let Ok(c) = Curve2::from_points(&points, self.tol, false) {
                Some(c)
            } else {
                None
            }
        }
    }

    pub fn normal_at(&self, l: f64) -> UnitVec2 {
        let r = self.at_length(l);

        if r.f <= self.tol {
            self.normal_at_vertex(r.i)
        } else if (r.f - 1.0).abs() <= self.tol {
            self.normal_at_vertex(r.i + 1)
        } else {
            self.line.edges()[r.i].normal.unwrap()
        }
    }

    pub fn direction_at(&self, l: f64) -> UnitVec2 {
        dir_from_normal(&self.normal_at(l))
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

        let mut lengths: Vec<f64> = vec![0.0];
        for e in line.edges().iter() {
            let d = dist(&line.points()[e.indices[0]], &line.points()[e.indices[1]]);
            lengths.push(d + lengths.last().unwrap_or(&0.0));
        }

        Ok(Curve2 {
            line,
            lengths,
            is_closed,
            tol,
        })
    }

    pub fn reversed(&self) -> Self {
        // TODO: this is ugly
        let r: Vec<Point2<f64>> = self.line.points().iter().rev().copied().collect_vec();
        Self::from_points(&r, self.tol, false).expect("Reverse failed but should not have")
    }

    pub fn edge_count(&self) -> usize {
        self.line.edges().len()
    }

    pub fn make_hull(&self) -> Option<ConvexPolygon<f64>> {
        ConvexPolygon::try_from_points(&self.line.points())
    }

    pub fn spanning_ray(&self, test_ray: &Ray<f64>) -> Option<SpanningRay> {
        spanning_ray(&self.line, test_ray)
    }

    pub fn max_ray_intersection(&self, ray: &Ray<f64>) -> Option<Point2<f64>> {
        max_intersection(&self.line, ray).map(|t| ray.point_at(t))
    }

    pub fn max_point_in_ray_direction(&self, ray: &Ray<f64>) -> Option<Point2<f64>> {
        let mut d = f64::MIN;
        let mut r = None;
        for p in self.line.points().iter() {
            let t = ray.projected_parameter(p);
            if t > d {
                d = t;
                r = Some(*p)
            }
        }
        r
    }
}

fn dir_from_normal(u: &UnitVec2) -> UnitVec2 {
    Isometry::rotation(std::f64::consts::PI * 0.5) * u
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::serialize::points_from_str;
    use approx::assert_relative_eq;
    use ncollide2d::na::{Point2, Vector2};
    use std::ffi::c_float;
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

    fn sample_points_scaled(p: &[(f64, f64)], f: f64) -> Vec<Point2<f64>> {
        p.iter().map(|(a, b)| Point2::new(*a * f, *b * f)).collect()
    }

    #[test]
    fn test_closest_point() {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();

        let p = curve.closest_point_and_edge(&Point2::new(2.0, 0.0));
        assert!(p.0 == 0 || p.0 == 1);
        assert_relative_eq!(1.0, p.1.x, epsilon = 1e-8);
        assert_relative_eq!(0.0, p.1.y, epsilon = 1e-8);
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
    #[test_case(-0.5, 0, 0.0)]
    #[test_case(0.0, 0, 0.0)]
    #[test_case(5.0, 3, 1.0)]
    #[test_case(2.0, 2, 0.0)]
    #[test_case(2.25, 2, 0.25)]
    fn test_lengths(l: f64, ei: usize, ef: f64) {
        let points = sample_points_scaled(&sample1(), 0.5);
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();

        let r = curve.at_length(l * 0.5);
        assert_eq!(ei, r.i);
        assert_relative_eq!(ef, r.f, epsilon = 1e-8);
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

    #[test_case(0.0, (-1.0, -1.0))]
    #[test_case(0.5, (0.0, -1.0))]
    #[test_case(1.0, (1.0, -1.0))]
    #[test_case(1.5, (1.0, 0.0))]
    #[test_case(2.0, (1.0, 1.0))]
    #[test_case(2.5, (0.0, 1.0))]
    #[test_case(3.0, (-1.0, 1.0))]
    #[test_case(3.5, (-1.0, 0.0))]
    #[test_case(4.0, (-1.0, -1.0))]
    fn test_normals_at_length_closed(l: f64, ec: (f64, f64)) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();
        let e = UnitVec2::new_normalize(Vector2::new(ec.0, ec.1));
        let n = curve.normal_at(l);

        assert_relative_eq!(e.x, n.x, epsilon = 1e-8);
        assert_relative_eq!(e.y, n.y, epsilon = 1e-8);
    }

    #[test_case(0.0, (0.0, -1.0))]
    #[test_case(0.5, (0.0, -1.0))]
    #[test_case(1.0, (1.0, -1.0))]
    #[test_case(1.5, (1.0, 0.0))]
    #[test_case(2.0, (1.0, 1.0))]
    #[test_case(2.5, (0.0, 1.0))]
    #[test_case(3.0, (0.0, 1.0))]
    fn test_normals_at_length_open(l: f64, ec: (f64, f64)) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, false).unwrap();
        let e = UnitVec2::new_normalize(Vector2::new(ec.0, ec.1));
        let n = curve.normal_at(l);

        assert_relative_eq!(e.x, n.x, epsilon = 1e-8);
        assert_relative_eq!(e.y, n.y, epsilon = 1e-8);
    }

    fn has_vertex(v: &Point2<f64>, c: &[Point2<f64>]) -> bool {
        for t in c.iter() {
            if dist(t, v) < 1e-6 {
                return true;
            }
        }
        false
    }

    #[test_case(0.0)]
    #[test_case(0.5)]
    #[test_case(0.75)]
    #[test_case(2.0)]
    #[test_case(2.1)]
    #[test_case(3.9)]
    fn test_distance_along(l: f64) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();
        let p = curve.point_at(l);
        let d = curve.distance_along(&p);

        assert_relative_eq!(l, d, epsilon = 1e-6);
    }

    #[test_case((0.1, 1.2), false, vec![1])] //             (0) |->  (1)  ->| (2)      (3)      O/C
    #[test_case((0.1, 2.2), false, vec![1, 2])] //          (0) |->  (1)  ->  (2)  ->| (3)      O/C
    #[test_case((0.7, 0.2), true, vec![1, 2, 3, 0])] //     (0)->||->(1)  ->  (2)  ->  (3)      C
    #[test_case((1.7, 1.2), true, vec![2, 3, 0, 1])] //     (0)  ->  (1)->||->(2)  ->  (3)      C
    #[test_case((2.7, 2.2), true, vec![3, 0, 1, 2])] //     (0)  ->  (1)  ->  (2)->||->(3) ->   C
    #[test_case((3.7, 3.2), true, vec![0, 1, 2, 3])] //     (0)  ->  (1)  ->  (2)  ->  (3)->||->C
    #[test_case((1.2, 0.7), true, vec![2, 3, 0])] //        (0)  ->| (1) |->  (2)  ->  (3) ->   C
    #[test_case((3.2, 0.7), true, vec![0])] //              (0)  ->| (1)      (2)      (3) ->|  C
    #[test_case((0.2, 3.7), true, vec![1, 2, 3])] //        (0) |->  (1)  ->  (2)  ->  (3) ->|  C
    #[test_case((0.1, 0.2), false, Vec::<usize>::new())] // (0) |->| (1)      (2)      (3)     O/C
    #[test_case((0.1, 0.2), true, Vec::<usize>::new())] //  (0) |->| (1)      (2)      (3)     O/C
    #[test_case((1.1, 1.8), false, Vec::<usize>::new())] // (0)      (1) |->| (2)      (3)     O/C
    #[test_case((1.1, 1.8), true, Vec::<usize>::new())] //  (0)      (1) |->| (2)      (3)     O/C
    #[test_case((3.1, 3.8), true, Vec::<usize>::new())] //  (0)      (1)      (2)      (3)|->| C
    fn test_portioning(l: (f64, f64), c: bool, i: Vec<usize>) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, c).unwrap();
        let p0 = curve.point_at(l.0);
        let p1 = curve.point_at(l.1);
        let result = curve.portion_between_lengths(l.0, l.1).unwrap();

        let e_l = if l.1 > l.0 {
            l.1 - l.0
        } else {
            curve.length() - (l.0 - l.1)
        };

        assert_relative_eq!(e_l, result.length(), epsilon = result.tol);
        assert_relative_eq!(p0.x, result.first_v().x, epsilon = result.tol);
        assert_relative_eq!(p0.y, result.first_v().y, epsilon = result.tol);
        assert_relative_eq!(p1.x, result.last_v().x, epsilon = result.tol);
        assert_relative_eq!(p1.y, result.last_v().y, epsilon = result.tol);

        for index in i {
            assert!(has_vertex(&points[index], result.line.points()));
        }
    }

    #[test]
    fn test_closest_point_search() {
        let points = points_from_str(include_str!("test_data/naca.txt"));
        let curve = Curve2::from_points(&points, 1e-4, false).unwrap();
        let query = Point2::new(-0.020415658216408356, 0.05938468508636321);

        let (edge, p) = curve.closest_point_and_edge(&query);
        assert_eq!(edge, 380);
        assert_relative_eq!(query.x, p.x, epsilon = 2e-6);
        assert_relative_eq!(query.y, p.y, epsilon = 2e-6);
    }

    #[test]
    fn test_naca_issue() {
        // l0: 0.061968421528170135
        // Lower:  0.06208009487909469, 0.0031649020793192567
        // Tested: 0.0595763602567424, 0.0020915926595290327
        // l1: 9.707314721439039
        // Upper:  -0.020415658216408356, 0.05938468508636321
        // Tested: -0.018854469810085348, 0.08048418445726035

        let points = points_from_str(include_str!("test_data/naca.txt"));
        let curve = Curve2::from_points(&points, 1e-4, false).unwrap();

        // let lower = Point2::new(0.06208009487909469, 0.0031649020793192567);
        let upper = Point2::new(-0.020415658216408356, 0.05938468508636321);

        let l = curve.distance_along(&upper);
        let p = curve.point_at(l);

        assert_relative_eq!(upper.x, p.x, epsilon = 2e-4);
        assert_relative_eq!(upper.y, p.y, epsilon = 2e-4);
    }
}
