use crate::geometry::shapes2::Circle2;
use ncollide2d::na::{Isometry2, Point2, Vector2};
use ncollide2d::shape::Ball;

use crate::airfoil::{Airfoil, CamberStation};
use crate::geometry::distances2::{dist, signed_angle};
use crate::geometry::line2::{intersect_rays, Line2};
use ncollide2d::query::Ray;
use std::error::Error;

const EPSILON: f64 = 1e-3;

fn deviation(p0: &Point2<f64>, p1: &Point2<f64>, test: &Point2<f64>) -> f64 {
    let ray = Ray::new(*p0, (p1 - p0).normalize());
    dist(&ray.projected_point(test), test)
}

/// An AirfoilGenerator is an entity which can generate the x, y position of the mean camber line
/// and the airfoil thickness at fractions of the chord. This provides the information necessary
/// for a generator to compute the airfoil surfaces.
pub trait AirfoilGenerator {
    /// Return a 2D point with the position of the camber line at a fraction from 0.0 to 1.0
    fn camber_line(&self, x: f64) -> Point2<f64>;

    /// Return the full thickness of the airfoil with respect to the camber line at a fraction from
    /// 0.0 to 1.0
    fn thickness(&self, x: f64) -> f64;

    fn station_at(&self, x: f64) -> CamberStation {
        let x0 = (x - EPSILON).max(0.0);
        let x1 = (x + EPSILON).min(1.0);

        let clx = self.camber_line(x);
        let cl0 = self.camber_line(x0);
        let cl1 = self.camber_line(x1);

        let d = Ray::new(clx, (cl1 - cl0).normalize());
        let n = d.turned();
        let t = self.thickness(x);

        CamberStation::new(clx, n.point_at(t / 2.0), n.point_at(-t / 2.0))
    }

    /// Adaptively generates airfoil information within a given tolerance
    fn generate(&self, tol: Option<f64>) -> Airfoil {
        let tol_value = tol.unwrap_or(1e-6);

        let mut stations: Vec<CamberStation> = Vec::new();
        let mut fractions: Vec<f64> = Vec::new();
        fractions.push(0.0);
        fractions.push(1.0);
        stations.push(self.station_at(0.0));
        stations.push(self.station_at(1.0));

        let mut index: usize = 0;
        while index < stations.len() - 1 {
            let x0 = &fractions[index];
            let x1 = &fractions[index + 1];

            let s0 = &stations[index];
            let s1 = &stations[index + 1];

            let x = (x0 + x1) / 2.0;
            let s = self.station_at(x);
            if deviation(&s0.upper, &s1.upper, &s.upper) < tol_value
                && deviation(&s0.lower, &s1.lower, &s.lower) < tol_value
                && deviation(&s0.camber, &s1.camber, &s.camber) < tol_value
            {
                index += 1;
            } else {
                fractions.insert(index + 1, x);
                stations.insert(index + 1, s);
            }
        }

        Airfoil::from_stations(&stations)
    }
}

/// A generator for a NACA 4-digit airfoil of the form MPTT, where M is the maximum camber P is the
/// location of the maximum camber, and TT is the maximum thickness of the airfoil as a fraction of
/// the chord.  For example, a NACA 2412 airfoil has a 2% camber at 40% of the chord and a max
/// thickness which is 12% of the chord length.
pub struct Naca4Digit {
    t: f64,
    chord_len: f64,
    m: f64,
    p: f64,
}

impl Naca4Digit {
    /// Create a new NACA 4 digit generator.
    ///
    /// # Arguments
    ///
    /// * `t_max` - the maximum thickness of the airfoil as a fraction of the chord length. For
    /// instance, on a NACA 2412 t_max should be 0.12
    ///
    /// * `chord_len` - the actual length of the airfoil chord
    ///
    /// * `max_camber` - The max camber as a fraction, for example on a NACA 2412 this value should
    /// be set to 0.02
    ///
    /// * `max_camber_chord` - The location of the max camber as a fraction of chord length. For
    /// example on a NACA 2412 this values should be 0.4
    pub fn new(t_max: f64, chord_len: f64, max_camber: f64, max_camber_chord: f64) -> Naca4Digit {
        Naca4Digit {
            t: t_max,
            chord_len,
            m: max_camber,
            p: max_camber_chord,
        }
    }
}

impl AirfoilGenerator for Naca4Digit {
    fn camber_line(&self, x: f64) -> Point2<f64> {
        let y = if self.p < 1e-6 {
            0.0
        } else if x < self.p {
            (self.m / self.p.powf(2.0)) * (2.0 * self.p * x - x.powf(2.0))
        } else {
            (self.m / (1.0 - self.p).powf(2.0))
                * ((1.0 - 2.0 * self.p) + 2.0 * self.p * x - x.powf(2.0))
        };

        Point2::new(x * self.chord_len, y * self.chord_len)
    }

    fn thickness(&self, x: f64) -> f64 {
        (2.0 * self.t * self.chord_len)
            * (1.485 * x.sqrt()
                + -0.630 * x
                + -1.758 * x.powf(2.0)
                + 1.4215 * x.powf(3.0)
                + -0.5075 * x.powf(4.0))
    }
}

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
        // The ends are parallel, so we should just do a full radius

        todo!()
    }
}

pub struct GeneratedAirfoil {
    pub mcl: Vec<Point2<f64>>,
    pub contour: Vec<Point2<f64>>,
    pub le: Circle2,
    pub te: Circle2,

    s0_indices: (usize, usize),
    s1_indices: (usize, usize),
    le_indices: (usize, usize),
    te_indices: (usize, usize),
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

        let s0_indices: (usize, usize) = (0, side0.len() - 1);
        side0.append(&mut tep);
        let te_indices: (usize, usize) = (s0_indices.1 + 1, side0.len() - 1);
        side0.append(&mut side1);
        let s1_indices: (usize, usize) = (te_indices.1 + 1, side0.len() - 1);
        side0.append(&mut lep);
        let le_indices: (usize, usize) = (s1_indices.1 + 1, side0.len() - 1);

        GeneratedAirfoil {
            mcl,
            contour: side0,
            le: lec,
            te: tec,
            s0_indices,
            s1_indices,
            le_indices,
            te_indices,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use serde::Deserialize;
    use test_case::test_case;

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

    #[test_case(1.000000, 0.001260)]
    #[test_case(0.840000, 0.021694)]
    #[test_case(0.680000, 0.038557)]
    #[test_case(0.520000, 0.051635)]
    #[test_case(0.360000, 0.059263)]
    #[test_case(0.200000, 0.057375)]
    #[test_case(0.040000, 0.032277)]
    fn test_naca_4_thickness(x: f64, e: f64) {
        let naca = Naca4Digit::new(0.12, 1.0, 0.0, 0.0);
        let result = naca.thickness(x);
        assert_relative_eq!(e * 2.0, result, epsilon = 1e-3);
    }

    #[test_case(1.000000, 0.001260)]
    #[test_case(0.840000, 0.021694)]
    #[test_case(0.680000, 0.038557)]
    #[test_case(0.520000, 0.051635)]
    #[test_case(0.360000, 0.059263)]
    #[test_case(0.200000, 0.057375)]
    #[test_case(0.040000, 0.032277)]
    fn test_naca_4_thickness_scaled(x: f64, e: f64) {
        let naca = Naca4Digit::new(0.12, 2.0, 0.0, 0.0);
        let result = naca.thickness(x);
        assert_relative_eq!(e * 4.0, result, epsilon = 1e-3);
    }

    #[test_case(1.0000, 0.0013)]
    #[test_case(0.9000, 0.0208)]
    #[test_case(0.7000, 0.0518)]
    #[test_case(0.5000, 0.0724)]
    #[test_case(0.3000, 0.0788)]
    #[test_case(0.2000, 0.0726)]
    #[test_case(0.1000, 0.0563)]
    fn test_naca_4_camber(x: f64, e: f64) {
        let naca = Naca4Digit::new(0.12, 1.0, 0.02, 0.4);
        let t = naca.thickness(x) / 2.0;
        let p = naca.camber_line(x);
        assert_relative_eq!(e, t + p.y, epsilon = 1e-3);
    }
}
