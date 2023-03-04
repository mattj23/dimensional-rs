use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ncollide2d::shape::Polyline;
use ncollide2d::query::Ray;
use ncollide2d::na::{ Point2, Vector2};

use dimensional_rs::geometry::polyline::{naive_ray_intersections, cleaned_polyline};
use dimensional_rs::airfoil::generate::{AirfoilGenerator, Naca4Digit};


fn polyline_intersections(line: &Polyline<f64>, rays: &[Ray<f64>]) -> Vec<f64> {
    let mut results = Vec::new();
    for ray in rays.iter() {
        results.append(&mut naive_ray_intersections(line, ray))
    }

    results
}

fn benchmark(c: &mut Criterion) {
    let naca = Naca4Digit::new(0.12, 10.0, 0.3, 0.4);
    let airfoil = naca.generate(Some(1e-4));
    let points = airfoil.to_outer_contour();
    let line = cleaned_polyline(&points, 1e-4);
    let ray_points: Vec<Point2<f64>> = (0..100).into_iter().map(|i| Point2::new(i as f64 / 10.0, 0.0)).collect();
    let rays: Vec<Ray<f64>> = ray_points.iter().map(|p| Ray::new(*p, Vector2::new(0.1, 0.9))).collect();

    c.bench_function("Polyline Intersections", |b| b.iter(|| polyline_intersections(black_box(&line), black_box(&rays))));
}


criterion_group!(benches, benchmark);
criterion_main!(benches);