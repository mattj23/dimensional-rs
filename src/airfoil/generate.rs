use crate::airfoil::{Airfoil, CamberStation};
use crate::geometry::distances2::deviation;
use crate::geometry::line2::Line2;
use ncollide2d::na::Point2;
use ncollide2d::query::Ray;

const EPSILON: f64 = 1e-3;

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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use test_case::test_case;

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
