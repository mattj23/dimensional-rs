use itertools::Itertools;
use ncollide2d::na::{Point2, Vector2};
use ncollide2d::query::Ray;
use ncollide2d::shape::Ball;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

#[derive(Serialize)]
#[serde(remote = "Point2<f64>")]
pub struct Point2f64 {
    x: f64,
    y: f64,
}

#[derive(Serialize)]
#[serde(remote = "Vector2<f64>")]
pub struct Vector2f64 {
    x: f64,
    y: f64,
}

#[derive(Serialize)]
#[serde(remote = "Ball<f64>")]
pub struct Ballf64 {
    radius: f64,
}

#[derive(Serialize)]
#[serde(remote = "Ray<f64>")]
pub struct Ray2f64 {
    #[serde(with = "Point2f64")]
    origin: Point2<f64>,

    #[serde(with = "Vector2f64")]
    dir: Vector2<f64>,
}

pub struct VectorList2f64 {
    vectors: Vec<Vector2<f64>>,
}

impl VectorList2f64 {
    pub fn from_vectors(vectors: &[Vector2<f64>]) -> Self {
        VectorList2f64 {
            vectors: vectors.to_vec(),
        }
    }

    pub fn from_points(points: &[Point2<f64>]) -> Self {
        VectorList2f64 {
            vectors: points.iter().map(|p| p.coords).collect(),
        }
    }
}

fn v2_to_str(v: &Vector2<f64>) -> String {
    format!("{},{}", v.x, v.y)
}

impl Serialize for VectorList2f64 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let s = format!("{}", self.vectors.iter().map(|v| v2_to_str(v)).format(";"));
        serializer.serialize_str(&s)
    }
}

pub fn points_to_string(points: &[Point2<f64>]) -> String {
    format!("{}", points.iter().map(|v| v2_to_str(&v.coords)).format(";"))
}

pub fn points_from_str(s: &str) -> Vec<Point2<f64>> {
    let mut raw = Vec::new();
    for token in s.split(";") {
        let mut pieces: [f64; 2] = [Default::default(); 2];
        for pair in token.split(",").take(2).enumerate() {
            pieces[pair.0] = pair
                .1
                .trim()
                .parse()
                .unwrap_or_else(|_| panic!("Couldn't parse a floating point value: '{}'", pair.1));
        }
        raw.push(Point2::new(pieces[0], pieces[1]));
    }
    raw
}

#[cfg(test)]
mod tests {
    use crate::serialize::VectorList2f64;
    use ncollide2d::na::Vector2;
    use serde_json;

    #[test]
    fn test_serialize_veclist() {
        let vecs = vec![
            Vector2::new(0.0, 1.0),
            Vector2::new(2.0, 3.0),
            Vector2::new(4.0, 5.0),
        ];
        let list = VectorList2f64::from_vectors(&vecs);

        let str = serde_json::to_string(&list).unwrap();
        assert_eq!("\"0,1;2,3;4,5\"", str);
    }
}
