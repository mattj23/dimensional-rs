use crate::geometry::shapes2::Circle2;
use ncollide2d::na::{Isometry2, Point2, Vector2};
use ncollide2d::shape::Ball;

use crate::geometry::distances2::{dist, signed_angle};
use crate::geometry::line2::{intersect_rays, Line2};
use ncollide2d::query::Ray;
use std::error::Error;

fn closest<'a>(i: &'a Point2<f64>, a: &'a Ray<f64>, b: &'a Ray<f64>) -> &'a Ray<f64> {
    if dist(i, &a.origin) < dist(i, &b.origin) {
        a
    } else {
        b
    }
}

fn edge_circle(
    a0: &Point2<f64>,
    a1: &Point2<f64>,
    b0: &Point2<f64>,
    b1: &Point2<f64>,
    n_points: usize,
) -> (Circle2, Vec<Point2<f64>>) {
    let a = Ray::new(a0.clone(), (a1 - a0).normalize());
    let b = Ray::new(b0.clone(), (b1 - b0).normalize());

    if let Some((a_t, _)) = intersect_rays(&a, &b) {
        // The two entities are not parallel
        let i0 = a.point_at(a_t);
        let ray = closest(&i0, &a, &b);

        // Get the middle ray
        let center_n = Isometry2::rotation(signed_angle(&a.dir, &b.dir) / 2.0) * a.dir;
        let center_ray = Ray::new(i0, center_n);
        let (cp0, _) = intersect_rays(&center_ray, &ray.turned()).unwrap();
        let center_point = center_ray.point_at(cp0);
        let radius = dist(&center_point, &ray.origin);
        let circle = Circle2::from_point(center_point, radius);

        let va = a.projected_point(&circle.center) - circle.center;
        let vb = b.projected_point(&circle.center) - circle.center;
        let angle = signed_angle(&va, &vb);

        let mut points = Vec::new();
        for i in 0..n_points + 1 {
            let pn: Vector2<f64> = Isometry2::rotation(angle * (i as f64) / (n_points as f64)) * va;
            let p: Point2<f64> = circle.center + pn;
            if dist(&p, &a0) > 1e-4 && dist(&p, &b0) > 1e-4 {
                points.push(p);
            }
        }

        (circle, points)
    } else {
        todo!()
    }
}

struct GeneratedAirfoil {
    mcl: Vec<Point2<f64>>,
    contour: Vec<Point2<f64>>,
    le: Circle2,
    te: Circle2,
}

impl GeneratedAirfoil {
    /// Generate an airfoil from a function which provides the mean camber line points and a
    /// another function which provides the thickness of the airfoil. Each function must take a
    /// value between 0.0 and 1.0 inclusive and return the point/thickness.  The parameter "n" will
    /// determine how many points will be sampled in the creation of the contour and mean camber
    /// line.
    ///
    /// The leading and trailing edge circles will be computed automatically from the tangencies
    /// left at the leading and trailing edge, which will extend beyond the 0.0 and 1.0 points of
    /// the mean camber line as they are retrieved from the functions.
    pub fn from_mcl_and_thickness<Fm, Ft>(
        f_mcl: Fm,
        f_thk: Ft,
        nm: usize,
        nc: usize,
    ) -> GeneratedAirfoil
    where
        Fm: Fn(f64) -> Point2<f64>,
        Ft: Fn(f64) -> f64,
    {
        let fractions: Vec<f64> = (0..nm + 1).map(|i| (i as f64) / (nm as f64)).collect();
        let mcl: Vec<Point2<f64>> = fractions.iter().map(|f| f_mcl(*f)).collect();
        let thk: Vec<f64> = fractions.iter().map(|f| f_thk(*f)).collect();

        let mut side0 = Vec::new();
        let mut side1 = Vec::new();

        for (i, value) in mcl.windows(2).enumerate() {
            let d = value[1] - value[0];
            let n = (Isometry2::rotation(std::f64::consts::PI / 2.0) * d).normalize();
            side0.push(value[0] + n * thk[i] / 2.0);
            side1.push(value[0] - n * thk[i] / 2.0);
        }

        // Get the very last thickness
        let last_d = mcl[nm] - mcl[nm - 1];
        let last_n = (Isometry2::rotation(std::f64::consts::PI / 2.0) * last_d).normalize();
        side0.push(mcl[nm] + last_n * thk[nm] / 2.0);
        side1.push(mcl[nm] - last_n * thk[nm] / 2.0);

        // Find the leading and trailing edge circles
        let (lec, mut lep) = edge_circle(&side0[0], &side0[1], &side1[0], &side1[1], nc);
        let (tec, mut tep) =
            edge_circle(&side1[nm], &side1[nm - 1], &side0[nm], &side0[nm - 1], nc);
        side1.reverse();
        lep.reverse();
        side0.append(&mut tep);
        side0.append(&mut side1);
        side0.append(&mut lep);

        GeneratedAirfoil {
            mcl,
            contour: side0,
            le: lec,
            te: tec,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;
    use approx::assert_relative_eq;

    #[derive(Deserialize)]
    struct TestGenAirfoilRaw {
        mcl: String,
        contour: String,
        le: String,
        te: String,
    }

    struct TestGenAirfoil {
        mcl: Vec<Point2<f64>>,
        contour: Vec<Point2<f64>>,
        le: Circle2,
        te: Circle2,
    }

    fn circle_from_string(s: &str) -> Result<(Point2<f64>, f64), Box<dyn Error>> {
        let mut values = Vec::new();
        for token in s.split(",") {
            values.push(token.parse()?);
        }

        Ok((Point2::new(values[0], values[1]), values[2]))
    }

    fn points_from_string(s: &str) -> Result<Vec<Point2<f64>>, Box<dyn Error>> {
        let mut raw_points: Vec<Point2<f64>> = Vec::new();
        for token in s.split(";") {
            let mut pieces: [f64; 2] = [Default::default(); 2];
            for pair in token.split(",").take(2).enumerate() {
                pieces[pair.0] = pair
                    .1
                    .parse()
                    .expect("Couldn't parse a floating point value");
            }
            raw_points.push(Point2::new(pieces[0], pieces[1]));
        }
        Ok(raw_points)
    }

    fn load_generated(raw_text: &str) -> Result<TestGenAirfoil, Box<dyn Error>> {
        let raw_data: TestGenAirfoilRaw = serde_json::from_str(&raw_text)?;
        let (le_point, le_radius) = circle_from_string(&raw_data.le)?;
        let (te_point, te_radius) = circle_from_string(&raw_data.te)?;
        Ok(TestGenAirfoil {
            mcl: points_from_string(&raw_data.mcl)?,
            contour: points_from_string(&raw_data.contour)?,
            le: Circle2::from_point(le_point, le_radius),
            te: Circle2::from_point(te_point, te_radius),
        })
    }

    fn tk(f: f64, t0: f64, t_max: f64, t1: f64, f_max: f64) -> f64 {
        if f < f_max {
            let f_ = f64::powf(1.0 - (f / f_max), 2.0);
            t0 * f_ + t_max * (1.0 - f_)
        } else {
            let f_ = f64::powf((f - f_max) / (1.0 - f_max), 2.0);
            t_max * (1.0 - f_) + t1 * f_
        }
    }

    #[test]
    fn test_generate() {
        let input = load_generated(include_str!("test_data/generated_airfoil_0.json"))
            .expect("Failed to deserialize");

        let circle = Circle2::new(0.0, 0.0, 3.0);
        let result = GeneratedAirfoil::from_mcl_and_thickness(
            |x| circle.point_at_angle(x * std::f64::consts::PI / 3.0),
            |x| tk(x, 0.25, 0.5, 0.05, 0.15),
            100,
            20,
        );

        for (i, v) in input.mcl.iter().enumerate() {
            assert_relative_eq!(v.x, result.mcl[i].x, epsilon = 1e-6);
            assert_relative_eq!(v.y, result.mcl[i].y, epsilon = 1e-6);
        }

        for (i, v) in input.contour.iter().enumerate() {
            let u = &result.contour[i];
            assert_relative_eq!(v.x, u.x, epsilon = 1e-6);
            assert_relative_eq!(v.y, u.y, epsilon = 1e-6);
        }


        assert_relative_eq!(input.le.center.x, result.le.center.x, epsilon = 1e-6);
        assert_relative_eq!(input.le.center.y, result.le.center.y, epsilon = 1e-6);
        assert_relative_eq!(input.le.ball.radius, result.le.ball.radius, epsilon = 1e-6);

        assert_relative_eq!(input.te.center.x, result.te.center.x, epsilon = 1e-6);
        assert_relative_eq!(input.te.center.y, result.te.center.y, epsilon = 1e-6);
        assert_relative_eq!(input.te.ball.radius, result.te.ball.radius, epsilon = 1e-6);
    }
}
