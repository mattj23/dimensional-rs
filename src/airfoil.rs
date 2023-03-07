use ncollide2d::na::Point2;

pub mod analyze;
mod common;
pub mod generate;

pub use analyze::analyze_airfoil;
pub use common::{AfParams, EdgeDetect};

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
    pub fn from_stations(stations: &[CamberStation]) -> Airfoil {
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
}

#[cfg(test)]
mod tests {}
