use dimensional_rs::airfoil::analyze::analyze_airfoil;
use dimensional_rs::airfoil::generate::{AirfoilGenerator, Naca4Digit};
use ncollide2d::na::Point2;
use std::fs::File;
use std::io::Write;

fn main() {
    let naca = Naca4Digit::new(0.12, 4.0, 0.3, 0.4);
    let airfoil = naca.generate(Some(1e-4));
    let outer_contour = airfoil.to_outer_contour();

    analyze_airfoil(&outer_contour, Some(1e-3));

    // write_points(&airfoil.upper, "data/side0.txt").expect("Failed writing file side0.txt");
    // write_points(&airfoil.lower, "data/side1.txt").expect("Failed writing file side1.txt");
    // write_points(&airfoil.to_outer_contour(), "data/contour.txt").expect("whoops");
}

// fn write_points(v: &Vec<Point2<f64>>, file_name: &str) -> std::io::Result<()> {
//     let mut file = File::create(file_name)?;
//     for p in v.iter() {
//         writeln!(file, "{}, {}", &p.x, &p.y).expect("Failed writing line");
//     }
//
//     Ok(())
// }
