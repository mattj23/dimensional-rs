use super::CamberStation;
use crate::closed_polyline::ClosedPolyline;
use crate::geometry::distances2::{
    closest_angle, deviation, dist, farthest_pair_indices, mid_point, signed_angle,
};
use crate::geometry::line2::Line2;
use crate::geometry::polyline::{cleaned_polyline, farthest_point_direction_distance, max_intersection, spanning_ray, SpanningRay};
use crate::geometry::shapes2::Circle2;
use ncollide2d::math::Isometry;
use ncollide2d::na::{Point2, Vector2};
use ncollide2d::query::{PointQuery, Ray};
use ncollide2d::shape::{ConvexPolygon, Polyline};
use std::fs::File;
use std::io::Write;
use crate::airfoil::analyze::EdgeDetect::Auto;

pub enum EdgeDetect {
    Auto,
    Circle,
    Square,
    Point
}

pub struct AfParams {
    pub tol: f64,
    pub leading: EdgeDetect,
    pub trailing: EdgeDetect
}

impl AfParams {
    pub fn new(tol: f64, leading: EdgeDetect, trailing: EdgeDetect) -> AfParams {
        AfParams {tol, leading, trailing}
    }

    pub fn default() -> AfParams{
        AfParams::new(1e-3, Auto, Auto)
    }
}

fn write_points(v: &[Point2<f64>], file_name: &str) -> std::io::Result<()> {
    let mut file = File::create(file_name)?;
    for p in v.iter() {
        writeln!(file, "{}, {}", &p.x, &p.y).expect("Failed writing line");
    }

    Ok(())
}
struct ExtractStation {
    spanning_ray: SpanningRay,
    upper: Point2<f64>,
    lower: Point2<f64>,
    circle: Circle2,
    thk: f64
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
            thk
        }
    }

    fn reversed(&self) -> ExtractStation {
        ExtractStation::new(self.spanning_ray.reversed(), self.lower.clone(), self.upper.clone(), self.circle.clone())
    }

    fn upper_dir(&self) -> Vector2<f64> {
        self.upper - self.circle.center
    }

    fn lower_dir(&self) -> Vector2<f64> {
        self.upper - self.circle.center
    }

    fn contact_angle(&self) -> f64 {
        closest_angle(&self.upper_dir(), &self.lower_dir())
    }

    /// Calculates the camber point between the upper and lower points
    fn cp(&self) -> Point2<f64> {
        mid_point(&self.upper, &self.lower)
    }

    fn contact_ray(&self) -> Ray<f64> {
        Ray::new(self.lower.clone(), self.upper - self.lower)
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

fn extract_station(line: &Polyline<f64>, ray: &SpanningRay, tol: Option<f64>) -> ExtractStation {
    let tol_value = tol.unwrap_or(1e-5);
    let mut positive = MeanSearchState::new(1.0, ray.at(1.0));
    let mut negative = MeanSearchState::new(0.0, ray.at(0.0));

    let mut working = ray.at(0.5);

    while (positive.fraction - negative.fraction) * ray.dir().norm() > tol_value {
        // println!("{}, {}", positive.fraction, negative.fraction);

        let fraction = (positive.fraction + negative.fraction) * 0.5;
        working = ray.at(fraction);

        let closest = line.project_point(&Isometry::identity(), &working, false);
        let to_closest = closest.point - working;
        let distance = dist(&working, &closest.point);
        if to_closest.dot(&ray.dir()) > 0.0 {
            positive.update(fraction, distance, closest.point.clone());
        } else {
            negative.update(fraction, distance, closest.point.clone());
        }
    }

    let circle = Circle2::from_point(
        ray.at((positive.fraction + negative.fraction) * 0.5),
        (positive.distance + negative.distance) * 0.5,
    );

    ExtractStation::new(ray.clone(), positive.point, negative.point, circle)
}

fn spanning_ray_advance(
    line: &Polyline<f64>,
    camber_ray: &Ray<f64>,
    d: f64,
    min_tol: f64,
) -> Option<SpanningRay> {
    let mut advance = 0.25;
    while advance >= 0.05 {
        let test_ray = Ray::new(camber_ray.point_at(d * advance), -camber_ray.orthogonal());
        if let Some(ray) = spanning_ray(&line, &test_ray) {
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
    line: &Polyline<f64>,
    s0: &ExtractStation,
    s1: &ExtractStation,
    outer_tol: f64,
    inner_tol: f64,
) -> Vec<ExtractStation> {
    let mut stations = Vec::new();

    let mid_ray = s0.spanning_ray.symmetry(&s1.spanning_ray);
    if let Some(ray) = spanning_ray(&line, &mid_ray) {
        let station = extract_station(&line, &ray, Some(inner_tol));

        if deviation(&s0.upper, &s1.upper, &station.upper) <= outer_tol
            && deviation(&s0.lower, &s1.lower, &station.lower) <= outer_tol
            && deviation(&s0.cp(), &s1.cp(), &station.cp()) <= outer_tol
        {
            stations.push(station);
        } else {
            let mut fwd = refine_between(&line, &s0, &station, outer_tol, inner_tol);
            let mut aft = refine_between(&line, &station, &s1, outer_tol, inner_tol);
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
    line: &Polyline<f64>,
    starting_ray: &SpanningRay,
    tol: Option<f64>,
) -> Vec<ExtractStation> {
    let outer_tol = tol.unwrap_or(1e-3);
    let inner_tol = outer_tol * 1e-1;

    let mut stations = Vec::new();
    let mut ray = starting_ray.clone();

    loop {
        let station = extract_station(&line, &ray, Some(inner_tol));

        if let Some(last_station) = stations.last() {
            let mut refined = refine_between(&line, last_station, &station, outer_tol, inner_tol);
            stations.append(&mut refined);
        }

        // Now we want to find the distance to the end of the airfoil in this direction. As we get
        // close to the leading or trailing edge, this distance will converge with the current
        // radius of the inscribed circle.
        let camber_ray = station.camber_ray();
        let to_farthest = farthest_point_direction_distance(&line, &camber_ray);
        let half_thickness = station.thk * 0.5;

        stations.push(station);

        if half_thickness * 1.1 >= to_farthest {
            break;
        }

        let advance_by = half_thickness.min(to_farthest);
        let advance_tol = outer_tol * 10.0;
        if let Some(next_ray) = spanning_ray_advance(&line, &camber_ray, advance_by, advance_tol) {
            ray = next_ray;
            continue;
        } else {
            println!("failed spanning ray");
            break;
        }
    }

    stations
}

pub fn analyze_airfoil(points: &Vec<Point2<f64>>, params: AfParams) {
    let line = cleaned_polyline(&points, params.tol);
    let hull = ConvexPolygon::try_from_points(&points).unwrap();

    write_points(&line.points(), "data/outer_contour.txt").unwrap();

    let (i0, i1) = farthest_pair_indices(&hull);
    let rough_chord = Ray::new( hull.points()[i0], hull.points()[i1] - hull.points()[i0] );

    // Create the initial intersections
    let mid_ray = Ray::new(rough_chord.point_at(0.5), rough_chord.orthogonal());
    let spanning = spanning_ray(&line, &mid_ray).expect("Failed on middle ray");

    let mut stations = extract_half_camber_line(&line, &spanning, None);
    let mut half = extract_half_camber_line(&line, &spanning.reversed(), None);
    stations.reverse();
    stations.append(&mut half.iter().skip(1).map(|s| s.reversed()).collect());

    let mut file = File::create("data/output.json").expect("Failed to create file");

    writeln!(file, "{{").unwrap();
    writeln!(file, "\"stations\": [").unwrap();

    for (i, s) in stations.iter().enumerate() {
        writeln!(file, "{{").unwrap();

        write_point("camber", &s.circle.center, &mut file);
        writeln!(file, ",").unwrap();

        write_point("upper", &s.upper, &mut file);
        writeln!(file, ",").unwrap();

        write_point("lower", &s.lower, &mut file);
        // writeln!(file, ",").unwrap();
        // writeln!(file, "\"rad\": {}", rads[i]).unwrap();

        writeln!(file, "}}").unwrap();

        if i < stations.len() - 1 {
            writeln!(file, ",").unwrap()
        }
    }
    writeln!(file, "]").unwrap();
    writeln!(file, "}}").unwrap();

    // write_points(&ints, "data/ints.txt").unwrap();
}

fn write_point(name: &str, point: &Point2<f64>, file: &mut File) {
    write!(file, "\"{}\": [{}, {}]", name, point.x, point.y).unwrap()
}