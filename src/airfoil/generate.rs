use crate::geometry::shapes2::Circle2;
use ncollide2d::na::Point2;
use ncollide2d::shape::Ball;

use std::error::Error;

struct GeneratedAirfoil {
    mcl: Vec<Point2<f64>>,
    contour: Vec<Point2<f64>>,
    le: Circle2,
    te: Circle2,
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;

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

    #[test]
    fn test_generate() {
        let input = load_generated(include_str!("test_data/generated_airfoil_0.json"))
            .expect("Failed to deserialize");

        assert_eq!(input.mcl, vec![Point2::<f64>::new(0.0, 0.1)]);
    }
}
