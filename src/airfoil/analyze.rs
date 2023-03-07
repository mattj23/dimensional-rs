use crate::geometry::curve2::Curve2;
use crate::geometry::distances2::{deviation, dist, farthest_pair_indices, mid_point};
use crate::geometry::line2::{intersect_rays, Line2};
use crate::geometry::polyline::{farthest_point_direction_distance, SpanningRay};
use crate::geometry::shapes2::Circle2;
use crate::serialize::{points_to_string, Point2f64, Vector2f64, VectorList2f64};
use ncollide2d::math::Isometry;
use ncollide2d::na::{Point2, Vector2};
use ncollide2d::query::{PointQuery, Ray};

use crate::airfoil::analyze::AdvanceResult::{CloseToEnd, FailedRay, NextRay};
use crate::airfoil::edges::EdgeDetect;
use serde::Serialize;
use serde_json;
use std::error::Error;
use std::fs::File;
use std::io::Write;

pub struct AfParams {
    pub tol: f64,
    pub leading: EdgeDetect,
    pub trailing: EdgeDetect,
}

impl Default for AfParams {
    fn default() -> Self {
        AfParams::new(1e-4, EdgeDetect::Auto, EdgeDetect::Auto)
    }
}

impl AfParams {
    pub fn auto_edges(tol: f64) -> Self {
        Self::new(tol, EdgeDetect::Auto, EdgeDetect::Auto)
    }

    pub fn new(tol: f64, leading: EdgeDetect, trailing: EdgeDetect) -> Self {
        AfParams {
            tol,
            leading,
            trailing,
        }
    }
}

#[derive(Serialize)]
struct InscribedCircle {
    spanning_ray: SpanningRay,

    #[serde(with = "Point2f64")]
    upper: Point2<f64>,

    #[serde(with = "Point2f64")]
    lower: Point2<f64>,
    circle: Circle2,
    thk: f64,
}

impl InscribedCircle {
    fn new(
        spanning_ray: SpanningRay,
        upper: Point2<f64>,
        lower: Point2<f64>,
        circle: Circle2,
    ) -> InscribedCircle {
        let thk = dist(&upper, &lower);
        InscribedCircle {
            spanning_ray,
            upper,
            lower,
            circle,
            thk,
        }
    }

    pub fn reversed(&self) -> InscribedCircle {
        InscribedCircle::new(
            self.spanning_ray.reversed(),
            self.lower,
            self.upper,
            self.circle.clone(),
        )
    }

    fn radius(&self) -> f64 {
        self.circle.ball.radius
    }

    fn center(&self) -> Point2<f64> {
        self.circle.center
    }

    /// Calculates the camber point between the upper and lower points
    fn cp(&self) -> Point2<f64> {
        mid_point(&self.upper, &self.lower)
    }

    fn contact_ray(&self) -> Ray<f64> {
        Ray::new(self.lower, self.upper - self.lower)
    }

    fn camber_ray(&self) -> Ray<f64> {
        Ray::new(self.center(), self.contact_ray().orthogonal().normalize())
    }

    fn interpolation_error(&self, s0: &Self, s1: &Self) -> f64 {
        deviation(&s0.upper, &s1.upper, &self.upper)
            .max(deviation(&s0.lower, &s1.lower, &self.lower))
            .max(deviation(&s0.center(), &s1.center(), &self.center()))
    }
}

struct MeanSearchState {
    fraction: f64,
    distance: f64,
    point: Point2<f64>,
}

impl MeanSearchState {
    fn new(fraction: f64, point: Point2<f64>) -> MeanSearchState {
        MeanSearchState {
            fraction,
            distance: 0.0,
            point,
        }
    }

    fn update(&mut self, fraction: f64, distance: f64, point: Point2<f64>) {
        self.distance = distance;
        self.fraction = fraction;
        self.point = point;
    }
}

/// Calculates the values of an inscribed circle station by using a maximum distance binary search
/// along a spanning ray through the curve. The search terminates when its motion falls below the
/// provided tolerance value.
fn calculate_inscribed(curve: &Curve2, ray: &SpanningRay, tol: f64) -> InscribedCircle {
    let mut positive = MeanSearchState::new(1.0, ray.at(1.0));
    let mut negative = MeanSearchState::new(0.0, ray.at(0.0));

    let mut working;

    while (positive.fraction - negative.fraction) * ray.dir().norm() > tol {
        let fraction = (positive.fraction + negative.fraction) * 0.5;
        working = ray.at(fraction);

        let closest = curve
            .line
            .project_point(&Isometry::identity(), &working, false);
        let to_closest = closest.point - working;
        let distance = dist(&working, &closest.point);
        if to_closest.dot(&ray.dir()) > 0.0 {
            positive.update(fraction, distance, closest.point);
        } else {
            negative.update(fraction, distance, closest.point);
        }
    }

    let circle = Circle2::from_point(
        ray.at((positive.fraction + negative.fraction) * 0.5),
        (positive.distance + negative.distance) * 0.5,
    );

    InscribedCircle::new(ray.clone(), positive.point, negative.point, circle)
}

enum AdvanceResult {
    NextRay(SpanningRay),
    CloseToEnd,
    FailedRay,
}

fn advance_ray(curve: &Curve2, station: &InscribedCircle) -> AdvanceResult {
    let camber_ray = station.camber_ray();
    let farthest = curve
        .max_point_in_ray_direction(&camber_ray)
        .expect("Failed to get farthest point while advancing ray");

    let d = dist(&farthest, &camber_ray.origin);

    if d - station.radius() < station.radius() * 0.5 {
        return CloseToEnd;
    }

    let mut f = 0.25;
    while f >= 0.05 {
        let next_center = camber_ray.point_at(f * station.radius());
        let test_ray = Ray::new(next_center, -camber_ray.orthogonal());
        if let Some(ray) = curve.spanning_ray(&test_ray) {
            if ray.dir().norm() < station.thk * 0.1 {
                f -= 0.05;
            } else {
                return NextRay(ray);
            }
        } else {
            f -= 0.05;
        }
    }

    FailedRay
}

fn refine_between(
    curve: &Curve2,
    s0: &InscribedCircle,
    s1: &InscribedCircle,
    outer_tol: f64,
    inner_tol: f64,
) -> Vec<InscribedCircle> {
    let mut stations = Vec::new();

    let mid_ray = s0.spanning_ray.symmetry(&s1.spanning_ray);
    if let Some(ray) = curve.spanning_ray(&mid_ray) {
        let station = calculate_inscribed(curve, &ray, inner_tol);

        if deviation(&s0.upper, &s1.upper, &station.upper) <= outer_tol
            && deviation(&s0.lower, &s1.lower, &station.lower) <= outer_tol
            && deviation(&s0.cp(), &s1.cp(), &station.cp()) <= outer_tol
        {
            stations.push(station);
        } else {
            let mut fwd = refine_between(curve, s0, &station, outer_tol, inner_tol);
            let mut aft = refine_between(curve, &station, s1, outer_tol, inner_tol);
            stations.append(&mut fwd);
            stations.push(station);
            stations.append(&mut aft);
        }
    }

    stations
}

fn refine_stations( curve: &Curve2, dest: &mut Vec<InscribedCircle>, stack: &mut Vec<InscribedCircle>, tol: f64) {
    while let Some(next) = stack.pop() {
        if let Some(last) = dest.last() {
            let test_ray = next.spanning_ray.symmetry(&last.spanning_ray);
            if let Some(ray) = curve.spanning_ray(&test_ray) {
                let mut mid = calculate_inscribed(curve, &ray, tol);
                // TODO: check the distance between the centers to make sure we're not stuck
                if mid.interpolation_error(&next, last) > tol {
                    // We are out of tolerance, we need to put next back on the stack and then put
                    // the mid ray on top of it and try again
                    stack.push(next);
                    stack.push(mid);
                } else {
                    // We are within tolerance, we can put the next station in the destination
                    dest.push(next);
                }
            }
        } else {
            dest.push(next);
        }
    }
}

/// Extracts the unambiguous portion of a mean camber line in the orthogonal direction to a
/// starting spanning ray. This function will terminate when it gets within a half thickness from
/// the furthest point in the camber line direction
fn extract_half_camber_line(
    curve: &Curve2,
    starting_ray: &SpanningRay,
    tol: Option<f64>,
) -> Vec<InscribedCircle> {
    let outer_tol = tol.unwrap_or(1e-3);
    let inner_tol = outer_tol * 1e-2;

    println!("Tolerances: {}, {}", outer_tol, inner_tol);

    let mut stations: Vec<InscribedCircle> = Vec::new();
    let mut refine_stack: Vec<InscribedCircle> = Vec::new();
    let mut ray = starting_ray.clone();

    loop {
        refine_stack.push(calculate_inscribed(curve, &ray, inner_tol));
        refine_stations(curve, &mut stations, &mut refine_stack, outer_tol);

        let station = &stations.last().expect("Station was not transferred");

        match advance_ray(curve, station) {
            NextRay(r) => {
                ray = r;
            }
            CloseToEnd => {
                break;
            }
            FailedRay => {
                println!("Failed to advance spanning ray");
                break;
            }
        };
    }

    stations
}

#[derive(Serialize)]
struct DebugData {
    stations: Vec<InscribedCircle>,
    line: VectorList2f64,
    // end0: VectorList2f64,
    // end1: VectorList2f64,
    #[serde(with = "Point2f64")]
    p0: Point2<f64>,
    #[serde(with = "Point2f64")]
    p1: Point2<f64>,
}

fn simple_end_intersection(curve: &Curve2, stations: &Vec<InscribedCircle>) -> Point2<f64> {
    let s1 = &stations[stations.len() - 1];
    let s0 = &stations[stations.len() - 2];
    let v = Ray::new(
        s1.circle.center,
        (s1.circle.center - s0.circle.center).normalize(),
    );
    // curve.max_ray_intersection(&v).expect("Failed on end intersection")
    curve
        .max_point_in_ray_direction(&v)
        .expect("Failed on end search")
}

pub fn analyze_airfoil(points: &[Point2<f64>], params: &AfParams) -> Result<(), Box<dyn Error>> {
    // Create the oriented curve (normals facing outward) and the first spanning ray
    let (curve, spanning) = create_curve_and_orient(points, params)?;

    // Build the unambiguous portions of the airfoil, away from the leading and trailing edges
    let mut stations = extract_half_camber_line(&curve, &spanning, Some(params.tol));
    let p0 = simple_end_intersection(&curve, &stations);

    let half = extract_half_camber_line(&curve, &spanning.reversed(), Some(params.tol));
    let p1 = simple_end_intersection(&curve, &half);

    stations.reverse();
    stations.append(&mut half.iter().skip(1).map(|s| s.reversed()).collect());

    // Now deal with the leading and trailing edges
    // let end0 = extract_edge_data(&curve, stations.first().unwrap(), false).unwrap();
    // let end1 = extract_edge_data(&curve, stations.last().unwrap(), true).unwrap();

    let mut file = File::create("data/output.json").expect("Failed to create file");
    let data = DebugData {
        stations,
        line: VectorList2f64::from_points(curve.line.points()),
        // end0: VectorList2f64::from_points(end0.line.points()),
        // end1: VectorList2f64::from_points(end1.line.points()),
        p0,
        p1,
    };
    let s = serde_json::to_string(&data).expect("Failed to serialize");
    write!(file, "{}", s)?;

    Ok(())
}

/// Given a full curve and the last station, portion out the
fn extract_edge_data(curve: &Curve2, station: &InscribedCircle, invert: bool) -> Option<Curve2> {
    let edge_dir = if invert {
        -station.camber_ray().dir.normalize()
    } else {
        station.camber_ray().dir.normalize()
    };

    let l0 = curve.distance_along(&station.lower);
    let l1 = curve.distance_along(&station.upper);

    // If the direction of the curve at l0 is in the same direction as the edge direction, we do
    // not need to swap the lengths while extracting the curve portion
    if curve.direction_at(l0).dot(&edge_dir) > 0.0 {
        curve.portion_between_lengths(l0, l1)
    } else {
        curve.portion_between_lengths(l1, l0)
    }
}

fn create_curve_and_orient(
    points: &[Point2<f64>],
    params: &AfParams,
) -> Result<(Curve2, SpanningRay), Box<dyn Error>> {
    let curve = Curve2::from_points(points, params.tol, false)?;
    let hull = curve.make_hull().expect("Failed to build convex hull");

    let (i0, i1) = farthest_pair_indices(&hull);
    let rough_chord = Ray::new(hull.points()[i0], hull.points()[i1] - hull.points()[i0]);

    // Create the initial intersection
    let mid_ray = Ray::new(rough_chord.point_at(0.5), rough_chord.orthogonal());

    let spanning = curve.spanning_ray(&mid_ray).expect("Failed on middle ray");
    let d = curve.distance_along(&spanning.origin());
    let n = curve.normal_at(d);

    if n.dot(&spanning.dir()) > 0.0 {
        Ok((curve.reversed(), spanning))
    } else {
        Ok((curve, spanning))
    }
}
