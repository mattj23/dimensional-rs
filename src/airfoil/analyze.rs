use super::CamberStation;
use crate::closed_polyline::ClosedPolyline;
use crate::geometry::distances2::{closest_angle, deviation, dist, farthest_pair_indices, mid_point, signed_angle};
use crate::geometry::line2::Line2;
use crate::geometry::polyline_intersections::{
    farthest_point_direction_distance, max_intersection, spanning_ray, SpanningRay,
};
use crate::geometry::shapes2::Circle2;
use ncollide2d::math::Isometry;
use ncollide2d::na::Point2;
use ncollide2d::query::{PointQuery, Ray};
use ncollide2d::shape::Polyline;
use std::fs::File;
use std::io::Write;

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
}

impl ExtractStation {
    fn new(
        spanning_ray: SpanningRay,
        upper: Point2<f64>,
        lower: Point2<f64>,
        circle: Circle2,
    ) -> ExtractStation {
        ExtractStation {
            spanning_ray,
            upper,
            lower,
            circle,
        }
    }

    fn contact_angle(&self) -> f64 {
       closest_angle(
            &(self.upper - self.circle.center),
            &(self.lower - self.circle.center),
        )
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
        // TODO: Check for when you're parallel to a square edge

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
) -> Option<SpanningRay> {
    let mut advance = 0.25;
    while advance >= 0.05 {
        let test_ray = Ray::new(camber_ray.point_at(d * advance), -camber_ray.orthogonal());
        if let Some(ray) = spanning_ray(&line, &test_ray) {
            return Some(ray);
        } else {
            advance -= 0.05;
        }
    }

    None
}

struct HalfExtractedAirfoil {
    stations: Vec<ExtractStation>,
    circle: Option<Circle2>,
    point: Option<Point2<f64>>,
}

impl HalfExtractedAirfoil {
    fn new(
        stations: Vec<ExtractStation>,
        circle: Option<Circle2>,
        point: Option<Point2<f64>>,
    ) -> HalfExtractedAirfoil {
        HalfExtractedAirfoil {
            stations,
            circle,
            point,
        }
    }
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

fn mcl_extract_rough(
    line: &Polyline<f64>,
    starting_ray: &SpanningRay,
    tol: Option<f64>,
) -> HalfExtractedAirfoil {
    /* This performs a rough extraction of the mean camber line on a closed contour starting at the
    starting_ray (which must be a spanning ray) and progressing in direction of the starting ray's
    orthogonal direction until reaching the end of the airfoil in the given direction.
     */

    let outer_tol = tol.unwrap_or(1e-3);
    let inner_tol = outer_tol * 1e-1;

    let mut stations = Vec::new();
    let mut ray = starting_ray.clone();

    loop {
        // Because we know that ray is a valid spanning ray, we can now find the center point along
        // the ray, which will yield the station.  However, in this case the station camber point
        // is the center of the inscribed circle, not the middle point between the upper and lower
        // surface points.
        let station = extract_station(&line, &ray, Some(inner_tol));

        // if let Some(last_station) = stations.last() {
        //     let mut refined = refine_between(&line, last_station, &station, outer_tol, inner_tol);
        //     stations.append(&mut refined);
        // }

        // The contact ray goes from the lower point to the upper point, such that its length is
        // the thickness and its center is on the mean camber line. The camber ray starts at the
        // center of the contact ray (on the mean camber line) and points in the direction of the
        // camber line at this station, with its length being the thickness of the airfoil.
        let contact_ray = station.contact_ray();
        let camber_ray = station.camber_ray();
        // println!("{:?}", camber_ray);
        let circle = station.circle.clone();

        // Now we want to find the distance to the end of the airfoil in this direction. As we get
        // close to the leading or trailing edge, this distance will converge with the current
        // radius of the inscribed circle.
        let to_farthest = farthest_point_direction_distance(&line, &camber_ray);
        let half_thickness = camber_ray.dir.norm() * 0.5;

        // If we are close to the leading or trailing edge on a closed section we will see two
        // things happening simultaneously: first, the half thickness will exceed the distance to
        // the farthest forward point, and second the intersection point of the camber ray with
        // the section itself will converge towards being on the inscribed circle of the station.

        // println!("{} >= {}", half_thickness, to_farthest);
        if half_thickness >= to_farthest {
            if let Some(t) = max_intersection(&line, &camber_ray) {
                // We have a forward intersection
                let end_point = camber_ray.point_at(t);

                /* The possible ending conditions are:
                   1. We have a contact circle with less than a 60 degree angle, indicating that
                       we probably have a leading edge radius
                   2. We have a point
                   3. We have a squared end
                */
                // println!("{:?}", circle.center);

                // Check for square ending condition first

                // Check for leading edge radius
                // println!("{}", station.contact_angle());
                if station.contact_angle().abs() < std::f64::consts::PI / 3.0 {
                    stations.push(station);
                    return HalfExtractedAirfoil::new(stations, Some(circle), Some(end_point));
                }

                // Check if we have converged
                if (dist(&end_point, &circle.center) - circle.ball.radius).abs() < outer_tol {}
            }
        }

        stations.push(station);

        let advance_by = half_thickness.min(to_farthest);
        if let Some(next_ray) = spanning_ray_advance(&line, &camber_ray, advance_by) {
            ray = next_ray;
            continue;
        } else {
            println!("failed spanning ray");
            break;
        }
    }

    HalfExtractedAirfoil::new(stations, None, None)
}

pub fn analyze_airfoil(points: &Vec<Point2<f64>>, join_tol: Option<f64>) {
    let contour = ClosedPolyline::new(&points, join_tol).expect("Failed to create polyline");

    write_points(&contour.line.points(), "data/outer_contour.txt").unwrap();

    let (i0, i1) = farthest_pair_indices(&contour.hull);
    let rough_chord = Ray::new(
        contour.hull.points()[i0],
        contour.hull.points()[i1] - contour.hull.points()[i0],
    );

    // Create the initial intersections
    let mid_ray = Ray::new(rough_chord.point_at(0.5), rough_chord.orthogonal());
    let mut spanning_ray = spanning_ray(&contour.line, &mid_ray).expect("Failed on middle ray");
    // spanning_ray = spanning_ray.reversed();

    let mut half = mcl_extract_rough(&contour.line, &spanning_ray, None);

    // // Debugging
    // let mut ints = Vec::new();
    //
    // // The
    // loop {
    //     ints.push(spanning_ray.point_at(0.0));
    //     ints.push(spanning_ray.point_at(1.0));
    //
    //     let station = find_center_point(&contour.line, &spanning_ray, Some(1e-5));
    //
    //     let contact_ray = Ray::new(station.lower.clone(), station.upper - station.lower);
    //     let motion_ray = Ray::new(contact_ray.point_at(0.5), contact_ray.orthogonal());
    //     let farthest = farthest_point_direction_distance(&contour.line, &motion_ray);
    //     // rads.push(farthest);
    //
    //     spanning_ray = contour
    //         .spanning_ray(&Ray::new(motion_ray.point_at(0.1), contact_ray.dir))
    //         .expect("Failed on calculating spanning ray");
    //     stations.push(station);
    //
    //     if farthest < contact_ray.dir.norm() {
    //         break;
    //     }
    //     // if counter > 10 {
    //     //     break;
    //     // }
    //     // counter+= 1;
    // }

    let mut file = File::create("data/output.json").expect("Failed to create file");

    writeln!(file, "{{").unwrap();
    writeln!(file, "\"stations\": [").unwrap();

    for (i, s) in half.stations.iter().enumerate() {
        writeln!(file, "{{").unwrap();

        write_point("camber", &s.circle.center, &mut file);
        writeln!(file, ",").unwrap();

        write_point("upper", &s.upper, &mut file);
        writeln!(file, ",").unwrap();

        write_point("lower", &s.lower, &mut file);
        // writeln!(file, ",").unwrap();
        // writeln!(file, "\"rad\": {}", rads[i]).unwrap();

        writeln!(file, "}}").unwrap();

        if i < half.stations.len() - 1 {
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
