use crate::geometry::curve2::Curve2;
use crate::geometry::distances2::{deviation, dist, farthest_pair_indices, mid_point};
use crate::geometry::line2::{intersect_rays, Line2};
use crate::geometry::polyline::{farthest_point_direction_distance, SpanningRay};
use crate::geometry::shapes2::Circle2;
use crate::serialize::{Point2f64, points_to_string, VectorList2f64};
use ncollide2d::math::Isometry;
use ncollide2d::na::Point2;
use ncollide2d::query::{PointQuery, Ray};

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
        AfParams::new(1e-3, EdgeDetect::Auto, EdgeDetect::Auto)
    }
}

impl AfParams {
    pub fn new(tol: f64, leading: EdgeDetect, trailing: EdgeDetect) -> Self {
        AfParams {
            tol,
            leading,
            trailing,
        }
    }
}

#[derive(Serialize)]
struct ExtractStation {
    spanning_ray: SpanningRay,

    #[serde(with = "Point2f64")]
    upper: Point2<f64>,

    #[serde(with = "Point2f64")]
    lower: Point2<f64>,
    circle: Circle2,
    thk: f64,
}

impl ExtractStation {
    fn new(
        spanning_ray: SpanningRay,
        upper: Point2<f64>,
        lower: Point2<f64>,
        circle: Circle2,
    ) -> ExtractStation {
        let thk = dist(&upper, &lower);
        ExtractStation {
            spanning_ray,
            upper,
            lower,
            circle,
            thk,
        }
    }

    pub fn reversed(&self) -> ExtractStation {
        ExtractStation::new(
            self.spanning_ray.reversed(),
            self.lower,
            self.upper,
            self.circle.clone(),
        )
    }

    // pub fn upper_dir(&self) -> Vector2<f64> {
    //     self.upper - self.circle.center
    // }

    // pub fn lower_dir(&self) -> Vector2<f64> {
    //     self.upper - self.circle.center
    // }

    // fn contact_angle(&self) -> f64 {
    //     closest_angle(&self.upper_dir(), &self.lower_dir())
    // }

    /// Calculates the camber point between the upper and lower points
    fn cp(&self) -> Point2<f64> {
        mid_point(&self.upper, &self.lower)
    }

    fn contact_ray(&self) -> Ray<f64> {
        Ray::new(self.lower, self.upper - self.lower)
    }

    fn camber_ray(&self) -> Ray<f64> {
        Ray::new(self.cp(), self.contact_ray().orthogonal().normalize())
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

fn extract_station(curve: &Curve2, ray: &SpanningRay, tol: Option<f64>) -> ExtractStation {
    let tol_value = tol.unwrap_or(1e-5);
    let mut positive = MeanSearchState::new(1.0, ray.at(1.0));
    let mut negative = MeanSearchState::new(0.0, ray.at(0.0));

    let mut working;

    while (positive.fraction - negative.fraction) * ray.dir().norm() > tol_value {
        // println!("{}, {}", positive.fraction, negative.fraction);

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

    ExtractStation::new(ray.clone(), positive.point, negative.point, circle)
}

fn spanning_ray_advance(
    curve: &Curve2,
    camber_ray: &Ray<f64>,
    d: f64,
    min_tol: f64,
) -> Option<SpanningRay> {
    let mut advance = 0.25;
    while advance >= 0.05 {
        let test_ray = Ray::new(camber_ray.point_at(d * advance), -camber_ray.orthogonal());
        if let Some(ray) = curve.spanning_ray(&test_ray) {
            if ray.dir().norm() < min_tol {
                advance -= 0.05;
            } else {
                return Some(ray);
            }
        } else {
            advance -= 0.05;
        }
    }

    None
}

fn refine_between(
    curve: &Curve2,
    s0: &ExtractStation,
    s1: &ExtractStation,
    outer_tol: f64,
    inner_tol: f64,
) -> Vec<ExtractStation> {
    let mut stations = Vec::new();

    let mid_ray = s0.spanning_ray.symmetry(&s1.spanning_ray);
    if let Some(ray) = curve.spanning_ray(&mid_ray) {
        let station = extract_station(curve, &ray, Some(inner_tol));

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

/// Extracts the unambiguous portion of a mean camber line in the orthogonal direction to a
/// starting spanning ray. This function will terminate when it gets within a half thickness from
/// the furthest point in the camber line direction
fn extract_half_camber_line(
    curve: &Curve2,
    starting_ray: &SpanningRay,
    tol: Option<f64>,
) -> Vec<ExtractStation> {
    let outer_tol = tol.unwrap_or(1e-3);
    let inner_tol = outer_tol * 1e-1;

    let mut stations: Vec<ExtractStation> = Vec::new();
    let mut ray = starting_ray.clone();

    loop {
        let station = extract_station(curve, &ray, Some(inner_tol));

        if let Some(last_station) = stations.last() {
            let d = station.cp() - last_station.cp();
            let mut refined = refine_between(curve, last_station, &station, outer_tol, inner_tol);
            stations.append(&mut refined);
        }

        // Now we want to find the distance to the end of the airfoil in this direction. As we get
        // close to the leading or trailing edge, this distance will converge with the current
        // radius of the inscribed circle.
        let camber_ray = station.camber_ray();
        let to_farthest = farthest_point_direction_distance(&curve.line, &camber_ray);
        let half_thickness = station.thk * 0.5;

        stations.push(station);

        if half_thickness * 1.1 >= to_farthest {
            break;
        }

        let (sr_dist, _) = intersect_rays(&camber_ray, &ray.ray())
            .expect("Failed measuring distance from camber ray to spanning ray");

        let advance_by = half_thickness.min(to_farthest) + sr_dist;
        let advance_tol = outer_tol * 10.0;
        if let Some(next_ray) = spanning_ray_advance(curve, &camber_ray, advance_by, advance_tol) {
            ray = next_ray;
            continue;
        } else {
            println!("failed spanning ray");
            break;
        }
    }

    stations
}

#[derive(Serialize)]
struct DebugData {
    stations: Vec<ExtractStation>,
    line: VectorList2f64,
    end0: VectorList2f64,
    end1: VectorList2f64,
}

pub fn analyze_airfoil(points: &[Point2<f64>], params: &AfParams) -> Result<(), Box<dyn Error>> {
    // Create the oriented curve (normals facing outward) and the first spanning ray
    let (curve, spanning) = create_curve_and_orient(points, params)?;

    // Build the unambiguous portions of the airfoil, away from the leading and trailing edges
    let mut stations = extract_half_camber_line(&curve, &spanning, None);
    let half = extract_half_camber_line(&curve, &spanning.reversed(), None);

    stations.reverse();
    stations.append(&mut half.iter().skip(1).map(|s| s.reversed()).collect());

    // Now deal with the leading and trailing edges
    // let end0 = extract_edge_data(&curve, stations.first().unwrap(), false).unwrap();
    // let end1 = extract_edge_data(&curve, stations.last().unwrap(), true).unwrap();

    let mut file = File::create("data/output.json").expect("Failed to create file");
    let data = DebugData {
        stations,
        line: VectorList2f64::from_points(curve.line.points()),
        end0: VectorList2f64::from_points(&Vec::new()),
        end1: VectorList2f64::from_points(&Vec::new()),

        // end0: VectorList2f64::from_points(end0.line.points()),
        // end1: VectorList2f64::from_points(end1.line.points()),
    };
    let s = serde_json::to_string(&data).expect("Failed to serialize");
    write!(file, "{}", s)?;

    Ok(())
}

/// Given a full curve and the last station, portion out the
fn extract_edge_data(curve: &Curve2, station: &ExtractStation, invert: bool) -> Option<Curve2> {
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
