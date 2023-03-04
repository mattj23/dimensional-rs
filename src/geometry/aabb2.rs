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
/// https://tavianator.com/2011/ray_box.html, except that I had to add branches to it because rust
/// seems to handle NAN propagation differently than C/C++
pub fn ray_intersects_aabb<N: RealField + Copy>(b: &AABB<N>, r: &Ray<N>, ix: N, iy: N) -> bool {
    let epsilon = N::default_epsilon() * N::from_f64(10.0).unwrap();
    // This is here because of issues getting the coincident slab edges to work
    if r.dir.x.abs() <= epsilon {
        r.origin.x + epsilon >= b.mins.x && r.origin.x - epsilon <= b.maxs.x
    } else if r.dir.y.abs() <= epsilon {
        r.origin.y + epsilon >= b.mins.y && r.origin.y - epsilon <= b.maxs.y
    } else {
        let tx1 = (b.mins.x - r.origin.x) * ix;
        let tx2 = (b.maxs.x - r.origin.x) * ix;
        let mut t_min = tx1.min(tx2);
        let mut t_max = tx1.max(tx2);

        let ty1 = (b.mins.y - r.origin.y) * iy;
        let ty2 = (b.maxs.y - r.origin.y) * iy;

        t_min = t_min.max(ty1.min(ty2));
        t_max = t_max.min(ty1.max(ty2));

        t_max >= t_min
    }
}

/// A visitor which traverses a BVH looking for intersections with a ray. Different from the
/// RayInterferenceCollector because it does not filter out intersections at negative parameters
pub struct RayVisitor<'a, N: 'a + RealField + Copy, T: 'a> {
    pub ray: &'a Ray<N>,
    pub collector: &'a mut Vec<T>,
    ix: N,
    iy: N,
}

impl<'a, N: RealField + Copy, T: Clone> RayVisitor<'a, N, T> {
    pub fn new(ray: &'a Ray<N>, buffer: &'a mut Vec<T>) -> RayVisitor<'a, N, T> {
        RayVisitor {
            ray,
            collector: buffer,
            ix: N::one() / ray.dir.x,
            iy: N::one() / ray.dir.y
        }
    }
}

impl<'a, N, T> Visitor<T, AABB<N>> for RayVisitor<'a, N, T>
where
    N: RealField + Copy,
    T: Clone,
{
    fn visit(&mut self, bv: &AABB<N>, t: Option<&T>) -> VisitStatus {
        if ray_intersects_aabb(bv, self.ray, self.ix, self.iy) {
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
    use ncollide2d::na::{Isometry2, Point2, Vector2};
    use test_case::test_case;

    #[test]
    fn test_ray_on_slab_0() {
        let aabb: AABB<f64> = AABB::new(Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));

        let ray: Ray<f64> = Ray::new(Point2::new(2.0, 1.0), Vector2::new(1.0, 0.0));
        let result = ray_intersects_aabb(&aabb, &ray, 1.0 / ray.dir.x, 1.0 / ray.dir.y);

        assert!(result);
    }

    #[test_case((2.0, 1.0, 1.0, 0.0), (0.0, 0.0))]
    #[test_case((2.0, 0.0, 1.0, 0.0), (0.0, 0.0))]
    fn test_ray_on_slab(a: (f64, f64, f64, f64), b: (f64, f64)) {
        let mut aabb: AABB<f64> = AABB::new(Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
        aabb.transform_by(&Isometry2::translation(b.0, b.1));

        let ray: Ray<f64> = Ray::new(Point2::new(a.0, a.1), Vector2::new(a.2, a.3));
        let result = ray_intersects_aabb(&aabb, &ray, 1.0 / ray.dir.x, 1.0 / ray.dir.y);

        assert!(result);
    }

    #[test]
    fn ray_aabb_test_intersect() {
        let aabb: AABB<f64> = AABB::new(Point2::new(0.0, 0.0), Point2::new(0.5, 0.5));
        let ray: Ray<f64> = Ray::new(Point2::new(0.25, 1.0), Vector2::new(0.0, 1.0));
        assert!(ray_intersects_aabb(&aabb, &ray, 1.0 / ray.dir.x, 1.0 / ray.dir.y));
    }

    #[test]
    fn ray_aabb_test_miss() {
        let aabb: AABB<f64> = AABB::new(Point2::new(0.0, 0.0), Point2::new(0.5, 0.5));
        let ray: Ray<f64> = Ray::new(Point2::new(1.25, 1.0), Vector2::new(0.0, 1.0));
        assert!(!ray_intersects_aabb(&aabb, &ray, 1.0 / ray.dir.x, 1.0 / ray.dir.y));
    }
}
