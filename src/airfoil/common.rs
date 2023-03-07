use crate::geometry::distances2::{deviation, dist, mid_point};
use crate::geometry::line2::Line2;
use crate::geometry::polyline::SpanningRay;
use crate::geometry::shapes2::Circle2;
use crate::serialize::Point2f64;
use ncollide2d::na::{Point2, Vector2};
use ncollide2d::query::Ray;
use serde::Serialize;
use std::mem::swap;

#[derive(Copy, Clone)]
pub enum EdgeDetect {
    Auto,
    Circle,
    Square,
    Point,
    Open,
}

/// This enum specifies the method for trying to detect the orientation of the leading edge on the
/// airfoil.
pub enum LeadingOrientation {
    /// This method is based on the principle that for most subsonic airfoils the point of max
    /// thickness is closer to the leading edge than to the trailing edge along the camber line
    TmaxForward,

    /// This method takes a direction vector and orients the leading edge such that the vector from
    /// the back to the front of the camber line has a positive dot product with that vector.
    Direction(Vector2<f64>),
}

pub struct AfParams {
    pub tol: f64,
    pub le_orientation: LeadingOrientation,
    pub leading: EdgeDetect,
    pub trailing: EdgeDetect,
}

impl Default for AfParams {
    fn default() -> Self {
        AfParams::new(
            1e-4,
            LeadingOrientation::TmaxForward,
            EdgeDetect::Auto,
            EdgeDetect::Auto,
        )
    }
}

impl AfParams {
    pub fn auto_edges(tol: f64) -> Self {
        Self::new(
            tol,
            LeadingOrientation::TmaxForward,
            EdgeDetect::Auto,
            EdgeDetect::Auto,
        )
    }

    pub fn new(
        tol: f64,
        le_orient: LeadingOrientation,
        leading: EdgeDetect,
        trailing: EdgeDetect,
    ) -> Self {
        AfParams {
            tol,
            le_orientation: le_orient,
            leading,
            trailing,
        }
    }
}

#[derive(Serialize)]
pub struct InscribedCircle {
    pub spanning_ray: SpanningRay,

    #[serde(with = "Point2f64")]
    pub upper: Point2<f64>,

    #[serde(with = "Point2f64")]
    pub lower: Point2<f64>,
    pub circle: Circle2,
    pub thk: f64,
}

impl InscribedCircle {
    pub fn new(
        spanning_ray: SpanningRay,
        upper: Point2<f64>,
        lower: Point2<f64>,
        circle: Circle2,
    ) -> InscribedCircle {
        let thk = dist(&upper, &lower);
        InscribedCircle {
            spanning_ray,
            upper,
            lower,
            circle,
            thk,
        }
    }

    pub fn reversed(&self) -> InscribedCircle {
        InscribedCircle::new(
            self.spanning_ray.reversed(),
            self.lower,
            self.upper,
            self.circle.clone(),
        )
    }

    pub fn reverse_in_place(&mut self) {
        self.spanning_ray = self.spanning_ray.reversed();
        swap(&mut self.upper, &mut self.lower)
    }

    pub fn radius(&self) -> f64 {
        self.circle.ball.radius
    }

    pub fn center(&self) -> Point2<f64> {
        self.circle.center
    }

    /// Calculates the camber point between the upper and lower points
    pub fn cp(&self) -> Point2<f64> {
        mid_point(&self.upper, &self.lower)
    }

    pub fn contact_ray(&self) -> Ray<f64> {
        Ray::new(self.lower, self.upper - self.lower)
    }

    pub fn camber_ray(&self) -> Ray<f64> {
        Ray::new(self.center(), self.contact_ray().orthogonal().normalize())
    }

    pub fn interpolation_error(&self, s0: &Self, s1: &Self) -> f64 {
        deviation(&s0.upper, &s1.upper, &self.upper)
            .max(deviation(&s0.lower, &s1.lower, &self.lower))
            .max(deviation(&s0.center(), &s1.center(), &self.center()))
    }
}
