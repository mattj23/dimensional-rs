use std::fs::File;
use std::io::Write;
use ncollide2d::na::Point2;
use dimensional_rs::airfoil::generate::{AirfoilGenerator, Naca4Digit};

fn main() {
    let naca = Naca4Digit::new(0.12, 4.0, 0.02, 0.4);
    let (a, b) = naca.generate(None);
    write_points(&a, "side0.txt").expect("Failed writing file side0.txt");
    write_points(&b, "side1.txt").expect("Failed writing file side1.txt");
}

fn write_points(v: &Vec<Point2<f64>>, file_name: &str) -> std::io::Result<()>  {
    let mut file = File::create(file_name)?;
    for p in v.iter() {
        writeln!(file, "{}, {}", &p.x, &p.y).expect("Failed writing line");
    }

    Ok(())
}

