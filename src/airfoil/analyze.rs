use super::{CamberStation};
use crate::closed_polyline::ClosedPolyline;
use crate::geometry::distances2::{dist, farthest_pair_indices};
use crate::geometry::line2::Line2;
use crate::geometry::polyline_intersections::{farthest_point_direction_distance, max_intersection, spanning_ray};
use ncollide2d::math::Isometry;
use ncollide2d::na::Point2;
use ncollide2d::query::{PointQuery, Ray};
use ncollide2d::shape::Polyline;
use std::fs::File;
use std::io::Write;
use crate::geometry::shapes2::Circle2;

fn write_points(v: &[Point2<f64>], file_name: &str) -> std::io::Result<()> {
    let mut file = File::create(file_name)?;
    for p in v.iter() {
        writeln!(file, "{}, {}", &p.x, &p.y).expect("Failed writing line");
    }

    Ok(())
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

fn find_center_point(line: &Polyline<f64>, ray: &Ray<f64>, tol: Option<f64>) -> CamberStation {
    let tol_value = tol.unwrap_or(1e-5);
    let mut positive = MeanSearchState::new(1.0, ray.point_at(1.0));
    let mut negative = MeanSearchState::new(0.0, ray.point_at(0.0));

    let mut working = ray.point_at(0.5);

    while (positive.fraction - negative.fraction) * ray.dir.norm() > tol_value {
        let fraction = (positive.fraction + negative.fraction) * 0.5;
        working = ray.point_at(fraction);

        let closest = line.project_point(&Isometry::identity(), &working, false);
        let to_closest = closest.point - working;
        let distance = dist(&working, &closest.point);
        if to_closest.dot(&ray.dir) > 0.0 {
            positive.update(fraction, distance, closest.point.clone());
        } else {
            negative.update(fraction, distance, closest.point.clone());
        }
    }

    CamberStation::new(working, positive.point, negative.point)
}


fn spanning_ray_advance(line: &Polyline<f64>, camber_ray: &Ray<f64>, d: f64) -> Option<Ray<f64>> {
    let mut advance = 0.25;
    let f = d / camber_ray.dir.norm();
    while advance >= 0.05 {
        let test_ray = Ray::new(camber_ray.point_at(f * advance), -camber_ray.orthogonal());
        if let Some(ray) = spanning_ray(&line, &test_ray) {
            return Some(ray);
        } else {
            advance -= 0.05;
        }
    };

    None
}

struct HalfExtractedAirfoil {
    stations: Vec<CamberStation>,
    circle: Option<Circle2>,
    point: Option<Point2<f64>>,
}

fn mcl_extract_rough(
    line: &Polyline<f64>,
    starting_ray: &Ray<f64>,
    tol: Option<f64>,
) -> HalfExtractedAirfoil {
    /* This performs a rough extraction of the mean camber line on a closed contour starting at the
    starting_ray (which must be a spanning ray) and progressing in direction of the starting ray's
    orthogonal direction until reaching the end of the airfoil in the given direction.
     */

    let outer_tol = tol.unwrap_or(1e-4);
    let inner_tol = outer_tol * 1e-1;

    let mut stations = Vec::new();
    let mut ray = starting_ray.clone();

    loop {
        // Because we know that ray is a valid spanning ray, we can now find the center point along
        // the ray, which will yield the station.  However, in this case the station camber point
        // is the center of the inscribed circle, not the middle point between the upper and lower
        // surface points.
        let station = find_center_point(&line, &ray, Some(inner_tol));

        // The contact ray goes from the lower point to the upper point, such that its length is
        // the thickness and its center is on the mean camber line. The camber ray starts at the
        // center of the contact ray (on the mean camber line) and points in the direction of the
        // camber line at this station, with its length being the thickness of the airfoil.
        let contact_ray = Ray::new(station.lower.clone(), station.upper - station.lower);
        let camber_ray = Ray::new(contact_ray.point_at(0.5), contact_ray.orthogonal());

        // Get the inscribed circle
        let circle = Circle2::from_point(station.camber.clone(), dist(&station.camber, &station.upper));
        stations.push(station);

        // Now we want to find the distance to the end of the airfoil in this direction. As we get
        // close to the leading or trailing edge, this distance will converge with the current
        // radius of the inscribed circle.
        let to_farthest = farthest_point_direction_distance(&line, &camber_ray);
        let half_thickness = camber_ray.dir.norm() * 0.5;

        // If we are close to the leading or trailing edge on a closed section we will see two
        // things happening simultaneously: first, the half thickness will exceed the distance to
        // the farthest forward point, and second the intersection point of the camber ray with
        // the section itself will converge towards being on the inscribed circle of the station.

        if half_thickness >= to_farthest {
            if let Some(t) = max_intersection(&line, &camber_ray) {
                // We have a forward intersection
                let end_point = camber_ray.point_at(t);

                // Check if we have converged
                if (dist(&end_point, &circle.center) - circle.ball.radius).abs() < outer_tol {
                    return HalfExtractedAirfoil{stations, circle: Some(circle), point: Some(end_point)};
                }
            }
        }

        // // circle center will drive
        // // trying to normally advance we
        // println!("{} >= {}", half_thickness, to_farthest);
        // if half_thickness >= to_farthest * 1.5 {
        //     println!("break");
        //     break
        // }

        let advance_by = half_thickness.max(to_farthest) * 0.5;
        if let Some(next_ray) = spanning_ray_advance(&line, &camber_ray, advance_by) {
            ray = next_ray;
            continue;
        } else {
            println!("failed spanning ray");
            break;
        }
    }

    HalfExtractedAirfoil{stations, circle: None, point: None}
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

        write_point("camber", &s.camber, &mut file);
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
