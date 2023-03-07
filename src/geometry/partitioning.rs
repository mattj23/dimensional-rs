use crate::geometry::aabb2::DistanceSearch;
use ncollide2d::na::{Point2, RealField};
use ncollide2d::partitioning::{VisitStatus, Visitor, BVH, BVT};
use std::collections::VecDeque;

#[derive(Debug)]
pub enum SearchType {
    Closest,
    Farthest,
}

pub struct PointVisitor<'a, N: 'a + RealField + Copy, T: 'a + Clone> {
    pub point: &'a Point2<N>,
    pub collector: &'a mut Vec<DistanceSearch<N, T>>,
    pub(crate) min_farthest: N,
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

pub trait BreadthFirst<T, BV> {
    fn bf_visit(&self, visitor: &mut impl Visitor<T, BV>);
}

impl<T, BV> BreadthFirst<T, BV> for BVT<T, BV> {
    fn bf_visit(&self, visitor: &mut impl Visitor<T, BV>) {
        let mut stack = VecDeque::new();

        if let Some(root) = self.root() {
            stack.push_back(root);

            while let Some(node) = stack.pop_front() {
                let content = self.content(node);

                match visitor.visit(content.0, content.1) {
                    VisitStatus::Continue => {
                        for i in 0..self.num_children(node) {
                            stack.push_back(self.child(i, node))
                        }
                    }
                    VisitStatus::ExitEarly => return,
                    VisitStatus::Stop => {}
                }
            }
        }
    }
}
