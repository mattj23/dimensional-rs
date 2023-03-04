use crate::geometry::distances2::dist;
use ncollide2d::bounding_volume::AABB;
use ncollide2d::na::{Point2, RealField};
use ncollide2d::partitioning::{VisitStatus, Visitor};
use ncollide2d::query::Ray;

pub fn aabb_closest_point<N: RealField + Copy>(b: &AABB<N>, p: &Point2<N>) -> Point2<N> {
    Point2::new(p.x.clamp(b.mins.x, b.maxs.x), p.y.clamp(b.mins.y, b.maxs.x))
}

pub fn aabb_farthest_point<N: RealField + Copy>(b: &AABB<N>, p: &Point2<N>) -> Point2<N> {
    let x = if (b.maxs.x - p.x).abs() > (b.mins.x - p.x).abs() {
        b.maxs.x
    } else {
        b.mins.x
    };

    let y = if (b.maxs.y - p.y).abs() > (b.mins.y - p.y).abs() {
        b.maxs.y
    } else {
        b.mins.y
    };

    Point2::new(x, y)
}

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
            iy: N::one() / ray.dir.y,
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

#[derive(Copy, Clone, Debug)]
pub struct DistanceSearch<N: RealField + Copy, T: Clone> {
    pub close: N,
    pub far: N,
    pub value: T,
}

impl<N: RealField + Copy, T: Clone> DistanceSearch<N, T> {
    fn new(close: N, far: N, value: T) -> DistanceSearch<N, T> {
        DistanceSearch { close, far, value }
    }
}

#[derive(Debug)]
pub enum SearchType {
    Closest,
    Farthest,
}

pub struct PointVisitor<'a, N: 'a + RealField + Copy, T: 'a + Clone> {
    pub point: &'a Point2<N>,
    pub collector: &'a mut Vec<DistanceSearch<N, T>>,
    min_farthest: N,
    max_closest: N,
    pub search_type: SearchType,
}

impl<'a, N: RealField + Copy, T: 'a + Clone> PointVisitor<'a, N, T> {
    pub fn new(
        point: &'a Point2<N>,
        buffer: &'a mut Vec<DistanceSearch<N, T>>,
        search_type: SearchType,
    ) -> PointVisitor<'a, N, T> {
        PointVisitor {
            point,
            collector: buffer,
            min_farthest: N::max_value().unwrap(),
            max_closest: N::min_value().unwrap(),
            search_type,
        }
    }
}

impl<'a, N, T> Visitor<T, AABB<N>> for PointVisitor<'a, N, T>
where
    N: RealField + Copy,
    T: Clone,
{
    fn visit(&mut self, bv: &AABB<N>, t: Option<&T>) -> VisitStatus {
        let close = dist(self.point, &aabb_closest_point(bv, self.point));
        let far = dist(self.point, &aabb_farthest_point(bv, self.point));

        match self.search_type {
            SearchType::Closest => {
                // If the closest distance to this bounding box is less than the smallest farthest
                // distance found so far, this is a candidate. If the farthest distance to this
                // bounding box is less than the smallest farthest distance found so far, then it
                // should be updated and any boxes which have a closest distance greater than the
                // new farthest distance should be ejected
                if far < self.min_farthest {
                    // Update and eject
                    self.min_farthest = far;
                    self.collector.retain(|d| d.close <= self.min_farthest);

                    // If we have a primitive we should collect it now
                    if let Some(t) = t {
                        self.collector
                            .push(DistanceSearch::new(close, far, t.clone()));
                    }

                    VisitStatus::Continue
                } else if close <= self.min_farthest {
                    // If we have a primitive we should collect it now
                    if let Some(t) = t {
                        self.collector
                            .push(DistanceSearch::new(close, far, t.clone()));
                    }

                    VisitStatus::Continue
                } else {
                    VisitStatus::Stop
                }
            }
            SearchType::Farthest => {
                todo!()
            }
        }

        // if ray_intersects_aabb(bv, self.ray, self.ix, self.iy) {
        //     if let Some(t) = t {
        //         self.collector.push(t.clone());
        //     }
        //
        //     VisitStatus::Continue
        // } else {
        //     VisitStatus::Stop
        // }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ncollide2d::na::{Isometry2, Point2, Vector2};
    use ncollide2d::partitioning::BVH;
    use ncollide2d::shape::Polyline;
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
        assert!(ray_intersects_aabb(
            &aabb,
            &ray,
            1.0 / ray.dir.x,
            1.0 / ray.dir.y
        ));
    }

    #[test]
    fn ray_aabb_test_miss() {
        let aabb: AABB<f64> = AABB::new(Point2::new(0.0, 0.0), Point2::new(0.5, 0.5));
        let ray: Ray<f64> = Ray::new(Point2::new(1.25, 1.0), Vector2::new(0.0, 1.0));
        assert!(!ray_intersects_aabb(
            &aabb,
            &ray,
            1.0 / ray.dir.x,
            1.0 / ray.dir.y
        ));
    }

    #[test]
    fn test_closest_point_visitor() {
        let poly = sample_polyline();
        let p: Point2<f64> = Point2::new(10.0, 0.0);
        let mut buffer = Vec::new();
        let mut v = PointVisitor::new(&p, &mut buffer, SearchType::Closest);

        poly.bvt().visit(&mut v);

        for i in buffer.iter() {
            //
        }

        assert!(false)
    }

    fn sample_polyline() -> Polyline<f64> {
        let points: Vec<Point2<f64>> = vec![
            Point2::new(5.0, 0.0),
            Point2::new(3.5, 0.9),
            Point2::new(4.0, 2.3),
            Point2::new(3.3, 3.3),
            Point2::new(2.7, 4.7),
            Point2::new(1.7, 6.4),
            Point2::new(0.0, 5.9),
            Point2::new(-1.5, 5.7),
            Point2::new(-3.7, 6.4),
            Point2::new(-5.3, 5.3),
            Point2::new(-6.4, 3.7),
            Point2::new(-7.1, 1.9),
            Point2::new(-7.3, 0.0),
            Point2::new(-7.8, -2.1),
            Point2::new(-6.3, -3.7),
            Point2::new(-5.7, -5.7),
            Point2::new(-3.7, -6.3),
            Point2::new(-1.7, -6.2),
            Point2::new(-0.0, -7.2),
            Point2::new(1.5, -5.6),
            Point2::new(2.4, -4.2),
            Point2::new(3.9, -3.9),
            Point2::new(4.9, -2.9),
            Point2::new(4.9, -1.3),
            Point2::new(5.0, 0.0),
        ];
        Polyline::new(points, None)
    }
}
