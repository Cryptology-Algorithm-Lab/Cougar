#![allow(warnings)]
mod BP;
mod Cougar;
mod Leopard;

// Test Codes
// use crate::BP::test_IPA::test_single;
// use crate::Leopard::test_IPA::test_single;
use crate::Cougar::test_IPA::test_single;

fn main() {
    // Enabling Single Thread for BP and Cougar
    use std::env;
    env::set_var("RAYON_NUM_THREADS", "1");

    // Recording the data
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create("your_file.txt").unwrap();

    for i in 10..21 {
        let (pt, vt, pc) = test_single(i);
        file.write_all((pt.to_string() + "\n").as_bytes()).unwrap();
        file.write_all((vt.to_string() + "\n").as_bytes()).unwrap();
        file.write_all((pc.to_string() + "\n").as_bytes()).unwrap();           
    }
    
}
