use dimensional_rs::airfoil::{analyze_airfoil, AfParams};
// use dimensional_rs::airfoil::generate::{AirfoilGenerator, Naca4Digit};
use dimensional_rs::serialize::points_from_str;
use std::fs;

fn main() {
    // let naca = Naca4Digit::new(0.12, 4.0, 0.3, 0.4);
    // let airfoil = naca.generate(Some(1e-4));
    // let outer_contour = airfoil.to_outer_contour();

    let data = fs::read_to_string("data/test.txt").unwrap();
    let points = points_from_str(&data);

    let result = analyze_airfoil(&points, &AfParams::auto_edges(1e-3));
}
