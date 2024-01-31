use std::fs::File;
use std::io::prelude::*;

pub fn main(){
    //file setup
    let mut render = File::create("image.ppm").expect("could not read ppm file");
    
    //img
    let image_width = 256;
    let image_height = 256;
    let starter = format!("P3\n{} {}\n225\n", image_width, image_height);

    write!(render, "{}", starter).expect("could not write to ppm");
    //rend
    for y in 0..image_height {
        for x in 0..image_width {
            let x_as_float = x as f32;
            let y_as_float = y as f32;
            let width_as_float = (image_width - 1) as f32;
            let height_as_float = (image_height - 1) as f32;
            let r = x_as_float / width_as_float;
            let g = y_as_float / height_as_float;

            let ir = (225.0 * r).round() as i32;
            let ig = (225.0 * g).round() as i32;

            let pixel = format!("{} {} 0\n", ir, ig); //rgb
            write!(render, "{}", pixel).expect("could not write to ppm");
        }
    }
}