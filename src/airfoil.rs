use ncollide2d::na::Vector2;

pub mod generate;

pub struct Airfoil {}

impl Airfoil {
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
