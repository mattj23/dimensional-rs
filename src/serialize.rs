use ncollide2d::na::{Point2, Vector2};
use ncollide2d::query::Ray;
use ncollide2d::shape::Ball;
use serde::Serialize;

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
