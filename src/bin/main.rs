use dimensional_rs::airfoil::analyze::{analyze_airfoil, AfParams};
use dimensional_rs::airfoil::generate::{AirfoilGenerator, Naca4Digit};

fn main() {
    let naca = Naca4Digit::new(0.12, 4.0, 0.3, 0.4);
    let airfoil = naca.generate(Some(1e-4));
    let outer_contour = airfoil.to_outer_contour();

    let result = analyze_airfoil(&outer_contour, &AfParams::default());
}

