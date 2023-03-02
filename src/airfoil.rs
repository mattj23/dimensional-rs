use ncollide2d::na::{Point2, Vector2};

pub mod analyze;
pub mod generate;

pub struct CamberStation {
    pub camber: Point2<f64>,
    pub upper: Point2<f64>,
    pub lower: Point2<f64>,
}

impl CamberStation {
    pub fn new(camber: Point2<f64>, upper: Point2<f64>, lower: Point2<f64>) -> CamberStation {
        CamberStation {
            camber,
            upper,
            lower,
        }
    }
}

pub struct Airfoil {
    pub camber: Vec<Point2<f64>>,
    pub upper: Vec<Point2<f64>>,
    pub lower: Vec<Point2<f64>>,
}

impl Airfoil {
    pub fn from_stations(stations: &Vec<CamberStation>) -> Airfoil {
        Airfoil {
            camber: stations.iter().map(|s| s.camber).collect(),
            upper: stations.iter().map(|s| s.upper).collect(),
            lower: stations.iter().map(|s| s.lower).collect(),
        }
    }

    pub fn to_outer_contour(&self) -> Vec<Point2<f64>> {
        let mut result = self.upper.to_vec();
        let mut lower = self.lower.to_vec();
        lower.reverse();
        result.append(&mut lower);
        result
    }

    pub fn analyze_from_points(points: &Vector2<f64>) -> Airfoil {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use generate::GeneratedAirfoil;

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
    fn test_airfoil_analyze() {}
}
