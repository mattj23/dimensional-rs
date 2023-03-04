use ncollide2d::bounding_volume::AABB;
use ncollide2d::na::{Point2, RealField};
use ncollide2d::partitioning::{VisitStatus, Visitor};
use ncollide2d::query::Ray;

/// Returns an array with the four corner points of an AABB
pub fn aabb_points<N: RealField + Copy>(b: &AABB<N>) -> [Point2<N>; 4] {
    [
        b.maxs,
        b.mins,
        Point2::new(b.mins.x, b.maxs.y),
        Point2::new(b.maxs.x, b.mins.y),
    ]
}

/// Detects whether or not an intersection exists between a ray and an aabb, regardless of the
/// direction or position of the intersection with respect to the ray origin.
///
/// Uses the fast, branchless method described by Tavian Barnes in his blog at
/// https://tavianator.com/2011/ray_box.html
pub fn ray_intersects_aabb<N: RealField + Copy>(b: &AABB<N>, r: &Ray<N>) -> bool {
    let x_inv = N::from_f64(1.0).unwrap() / r.dir.x;
    let y_inv = N::from_f64(1.0).unwrap() / r.dir.y;
    let tx1 = (b.mins.x - r.origin.x) * x_inv;
    let tx2 = (b.maxs.x - r.origin.x) * x_inv;
    let mut t_min = tx1.min(tx2);
    let mut t_max = tx1.max(tx2);

    let ty1 = (b.mins.y - r.origin.y) * y_inv;
    let ty2 = (b.maxs.y - r.origin.y) * y_inv;

    t_min = t_min.max(ty1.min(ty2));
    t_max = t_max.min(ty1.max(ty2));

    t_max >= t_min
}

/// A visitor which traverses a BVH looking for intersections with a ray. Different from the
/// RayInterferenceCollector because it does not filter out intersections at negative parameters
pub struct RayVisitor<'a, N: 'a + RealField + Copy, T: 'a> {
    pub ray: &'a Ray<N>,
    pub collector: &'a mut Vec<T>,
}

impl<'a, N: RealField + Copy, T: Clone> RayVisitor<'a, N, T> {
    pub fn new(ray: &'a Ray<N>, buffer: &'a mut Vec<T>) -> RayVisitor<'a, N, T> {
        RayVisitor {
            ray,
            collector: buffer,
        }
    }
}

impl<'a, N, T> Visitor<T, AABB<N>> for RayVisitor<'a, N, T>
where
    N: RealField + Copy,
    T: Clone,
{
    fn visit(&mut self, bv: &AABB<N>, t: Option<&T>) -> VisitStatus {
        if ray_intersects_aabb(bv, self.ray) {
            if let Some(t) = t {
                self.collector.push(t.clone());
            }

            VisitStatus::Continue
        } else {
            VisitStatus::Stop
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ncollide2d::na::{Point2, Vector2};

    #[test]
    fn ray_aabb_test_intersect() {
        let aabb: AABB<f64> = AABB::new(Point2::new(0.0, 0.0), Point2::new(0.5, 0.5));
        let ray: Ray<f64> = Ray::new(Point2::new(0.25, 1.0), Vector2::new(0.0, 1.0));
        assert!(ray_intersects_aabb(&aabb, &ray));
    }

    #[test]
    fn ray_aabb_test_miss() {
        let aabb: AABB<f64> = AABB::new(Point2::new(0.0, 0.0), Point2::new(0.5, 0.5));
        let ray: Ray<f64> = Ray::new(Point2::new(1.25, 1.0), Vector2::new(0.0, 1.0));
        assert!(!ray_intersects_aabb(&aabb, &ray));
    }
}
