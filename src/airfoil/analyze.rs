use crate::geometry::curve2::Curve2;
use crate::geometry::distances2::{deviation, dist, farthest_pair_indices};
use crate::geometry::line2::Line2;
use crate::geometry::polyline::SpanningRay;
use crate::geometry::shapes2::Circle2;
use crate::serialize::VectorList2f64;
use ncollide2d::math::Isometry;
use ncollide2d::na::{Point2, Vector2};
use ncollide2d::query::{PointQuery, Ray};

use super::common::{AfParams, EdgeDetect, InscribedCircle};
use crate::errors::InvalidGeometry;
use serde::Serialize;
use serde_json;
use std::error::Error;
use std::fs::File;
use std::io::Write;

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
        return AdvanceResult::CloseToEnd;
    }

    let mut f = 0.25;
    while f >= 0.05 {
        let next_center = camber_ray.point_at(f * station.radius());
        let test_ray = Ray::new(next_center, -camber_ray.orthogonal());
        if let Some(ray) = curve.spanning_ray(&test_ray) {
            if ray.dir().norm() < station.thk * 0.1 {
                f -= 0.05;
            } else {
                return AdvanceResult::NextRay(ray);
            }
        } else {
            f -= 0.05;
        }
    }

    AdvanceResult::FailedRay
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

fn refine_stations(
    curve: &Curve2,
    dest: &mut Vec<InscribedCircle>,
    stack: &mut Vec<InscribedCircle>,
    tol: f64,
) {
    while let Some(next) = stack.pop() {
        if let Some(last) = dest.last() {
            let test_ray = next.spanning_ray.symmetry(&last.spanning_ray);
            if let Some(ray) = curve.spanning_ray(&test_ray) {
                let mid = calculate_inscribed(curve, &ray, tol);
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
/// starting spanning ray. This function will terminate when it gets close to the farthest point in
/// the camber line direction.
fn extract_half_camber_line(
    curve: &Curve2,
    starting_ray: &SpanningRay,
    tol: Option<f64>,
) -> Vec<InscribedCircle> {
    let outer_tol = tol.unwrap_or(1e-3);
    let inner_tol = outer_tol * 1e-2;

    let mut stations: Vec<InscribedCircle> = Vec::new();
    let mut refine_stack: Vec<InscribedCircle> = Vec::new();
    let mut ray = starting_ray.clone();

    loop {
        refine_stack.push(calculate_inscribed(curve, &ray, inner_tol));
        refine_stations(curve, &mut stations, &mut refine_stack, outer_tol);

        let station = &stations.last().expect("Station was not transferred");

        match advance_ray(curve, station) {
            AdvanceResult::NextRay(r) => {
                ray = r;
            }
            AdvanceResult::CloseToEnd => {
                break;
            }
            AdvanceResult::FailedRay => {
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
    // #[serde(with = "Point2f64")]
    // p0: Point2<f64>,
    // #[serde(with = "Point2f64")]
    // p1: Point2<f64>,
}

fn reverse_stations(stations: &mut [InscribedCircle]) {
    stations.iter_mut().for_each(|i| i.reverse_in_place());
}

pub fn analyze_airfoil(points: &[Point2<f64>], params: &AfParams) -> Result<(), Box<dyn Error>> {
    // Create the oriented curve (normals facing outward) and the first spanning ray
    let (curve, spanning) = create_curve_and_orient(points, params)?;

    // Build the unambiguous portions of the airfoil, away from the leading and trailing edges
    let mut stations = extract_half_camber_line(&curve, &spanning, Some(params.tol));
    let half = extract_half_camber_line(&curve, &spanning.reversed(), Some(params.tol));
    stations.reverse();
    stations.append(&mut half.iter().skip(1).map(|s| s.reversed()).collect());

    if stations.is_empty() {
        return Err(InvalidGeometry::GeometricOpFailed.into());
    }

    // Attempt to detect the orientation of the leading edge

    // Now deal with the leading and trailing edges
    if let Some(end0) = extract_edge_data(&curve, stations.first().unwrap(), false) {
        // compute_edge(&curve, stations.first().unwrap(), params.)
    }
    // let end0 = extract_edge_data(&curve, stations.first().unwrap(), false).unwrap();
    // let end1 = extract_edge_data(&curve, stations.last().unwrap(), true).unwrap();

    let mut file = File::create("data/output.json").expect("Failed to create file");
    let data = DebugData {
        stations,
        line: VectorList2f64::from_points(curve.line.points()),
        // end0: VectorList2f64::from_points(end0.line.points()),
        // end1: VectorList2f64::from_points(end1.line.points()),
    };
    let s = serde_json::to_string(&data).expect("Failed to serialize");
    write!(file, "{}", s)?;

    Ok(())
}

fn le_is_oriented_by_tmax(curve: &Curve2) -> bool {
    todo!()
}

fn le_is_oriented_by_direction(curve: &Curve2, dir: &Vector2<f64>) -> bool {
    todo!()
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

fn compute_edge(curve: &Curve2, station: &InscribedCircle, edge_type: &EdgeDetect, tol: f64) {}
