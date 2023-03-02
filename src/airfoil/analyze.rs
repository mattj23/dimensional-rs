use super::{Airfoil, CamberStation};
use crate::closed_polyline::ClosedPolyline;
use crate::geometry::distances2::{dist, farthest_pair_indices};
use crate::geometry::line2::Line2;
use ncollide2d::na::Point2;
use ncollide2d::query::{PointQuery, Ray};
use ncollide2d::shape::Polyline;
use std::fs::File;
use std::io::Write;
use ncollide2d::math::Isometry;

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

    fn update(&mut self, fraction:f64, distance: f64, point: Point2<f64>) {
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

fn farthest_point_direction(line: &Polyline<f64>, ray: &Ray<f64>) -> f64 {
    let mut farthest: f64 = 0.0;
    for v in line.points().iter() {
        let d = v - ray.origin;
        if ray.dir.dot(&d) > 0.0 {
            farthest = farthest.max(d.norm());
        }
    }

    farthest
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
    let mut spanning_ray = contour.spanning_ray(&mid_ray).expect("Failed on middle ray");
    spanning_ray = spanning_ray.transform_by(&Isometry::rotation(std::f64::consts::PI));

    let mut stations = Vec::new();

    // Debugging
    let mut ints = Vec::new();
    // let mut rads = Vec::new();

    let mut counter = 0;
    loop {
        ints.push(spanning_ray.point_at(0.0));
        ints.push(spanning_ray.point_at(1.0));

        let station = find_center_point(&contour.line, &spanning_ray, Some(1e-5));

        let contact_ray = Ray::new(station.lower.clone(), station.upper - station.lower);
        let motion_ray = Ray::new(contact_ray.point_at(0.5), contact_ray.orthogonal());
        let farthest = farthest_point_direction(&contour.line, &motion_ray);
        // rads.push(farthest);

        spanning_ray = contour.spanning_ray(&Ray::new(motion_ray.point_at(0.5), contact_ray.dir))
            .expect("Failed on calculating spanning ray");
        stations.push(station);

        if farthest < contact_ray.dir.norm() * 0.5 {
            break;
        }
        if counter > 1 {
            break;
        }
        counter+= 1;
    }


    let mut file = File::create("data/output.json").expect("Failed to create file");

    writeln!(file, "{{").unwrap();
    writeln!(file, "\"stations\": [").unwrap();

    for (i, s) in stations.iter().enumerate() {
        writeln!(file, "{{").unwrap();

        write_point("camber", &s.camber, &mut file);
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

    write_points(&ints, "data/ints.txt").unwrap();
}


fn write_point(name: &str, point: &Point2<f64>, file: &mut File)  {
    write!(file, "\"{}\": [{}, {}]", name, point.x, point.y).unwrap()
}