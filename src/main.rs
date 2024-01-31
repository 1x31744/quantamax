extern crate sdl2;

use cond_utils::Between;
use sdl2::pixels::Color;
use std::fs::File;
use std::io::prelude::*;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::render::Canvas;
use sdl2::video::Window;
use core::panic;
use std::time::Duration;
use sdl2::rect::Point;

const WIDTH: u32 = 1000;
const HEIGHT: u32 = 500;

#[derive(PartialEq)]
struct Sphere {
    pub radius: u32,
    pub center: (i32, i32, i32),
    pub color: (u8, u8, u8),

}
impl Sphere {
    pub fn new(radius: u32, center: (i32, i32, i32), color: (u8, u8, u8)) -> Sphere {
        Sphere {
            radius: radius,
            center: center,
            color: color,
        }
    }
}


fn put_pixel(canvas: &mut Canvas<Window>,x: i32, y:i32, color: (u8, u8, u8)) {
    let window_width_half: i32 = (WIDTH/2) as i32;
    let window_height_half: i32 = (HEIGHT/2) as i32;
    //error handling for draw position out of range
    if x > window_width_half || x < -window_width_half {
        println!("invalid x = {}", x);
        panic!("x co-ordinate outside range")
    }
    if y > window_height_half || y < -window_height_half {
        println!("invalid y = {}", y);
        panic!("y co-ordinate outside range")
    }
    //convert cartesian to computer
    let canvas_x: i32 = (x + window_width_half) as i32;
    let canvas_y: i32 = (window_height_half - y) as i32;
    //draw pixel
    canvas.set_draw_color(Color::RGB(color.0, color.1, color.2));
    let point = Point::new(canvas_x, canvas_y);
    canvas.draw_point(point).expect("could not draw point");

    




}

// TODO: replace all tuples with lists as tuples are bad practice for same typing
// TODO: finish dot product function and then finish sphere ray intersection function.

pub fn main() {
    let mut debug = File::create("debug.txt").expect("could not read debug file");
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("quantamax", WIDTH, HEIGHT)
        .position_centered()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();
    let mut event_pump = sdl_context.event_pump().unwrap();

    let camera_origin = (0,0,0);
    let window_width_half: i32 = (WIDTH/2) as i32;
    let window_height_half: i32 = (HEIGHT/2) as i32;
    let mut sphere = Sphere::new(5, (0,0, 50), (225,0,0));
    'running: loop {
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..} |
                Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                    break 'running
                },
                _ => {}
            }
        }
        // The rest of the game loop goes here...

        canvas.set_draw_color(Color::RGB(225, 225, 225));
        canvas.clear();
        /* 
        let (render_width, render_height) = canvas.output_size().unwrap();
        let mut x_as_float: f32;
        let mut y_as_float: f32;
        let mut width_as_float: f32;
        let mut height_as_float: f32;
        let mut r: f32;
        let mut g: f32;
        let mut ir: u8;
        let mut ig: u8;
        
        
        for y in 0..render_height {
            for x in 0..render_width {
                x_as_float = x as f32;
                y_as_float = y as f32;
                width_as_float = (render_width - 1) as f32;
                height_as_float = (render_height - 1) as f32;
                r = x_as_float / width_as_float;
                g = y_as_float / height_as_float;

                ir = (225.0 * r).round() as u8;
                ig = (225.0 * g).round() as u8;
                //println!("{} {}",ir, ig);
                canvas.set_draw_color(Color::RGB(ir, ig, 0));
                let point = Point::new(x.try_into().unwrap(), y.try_into().unwrap());
                canvas.draw_point(point).expect("could not draw point");
            }
        }
        */

        for x in -window_width_half..window_width_half {
            for y in -window_height_half..window_height_half {
                let view = canvas_to_viewport(x, y, WIDTH as i32, HEIGHT as i32, 1);
                //let debug_text = format!("x, y = {}, {} \n view = {:?} \n\n", x, y, viewport);
                //write!(debug, "{}", debug_text);
                let color = trace_ray(camera_origin, view, 1, 100, &mut sphere);
           }
           //write!(debug, "{}", "new line");
        }
        println!("done");
        put_pixel(&mut canvas, 1, 1, (225,0,0));
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60));
    }
}

fn trace_ray(camera_origin: (i32, i32, i32), view: (i32, i32, i32), tmin: i32, tmax: i32, sphere: &mut Sphere) -> (u8, u8, u8) {
    let mut closest_sphere: Option<&mut Sphere> = None;
    let mut closest_t: i32 = 100000000;
    //for sphere in spheres
    let intersects : (i32, i32) = intersect_ray_sphere(camera_origin, view, sphere); //two intersections
    if intersects.0.within(tmin, tmax) && intersects.0 < closest_t {
        closest_sphere = Some(sphere)
    }
    if intersects.1.within(tmin, tmax) && intersects.1 < closest_t {
        closest_sphere = Some(sphere)
    }
    if closest_sphere == None {
        return (0,0,0)
    }
    return closest_sphere.unwrap().color;
}

fn intersect_ray_sphere(camera_origin: (i32, i32, i32), view: (i32, i32, i32), sphere: &mut Sphere) -> (i32, i32) {
    let radius = sphere.radius;
    let camera_to_center = (camera_origin.0 - sphere.center.0, camera_origin.1 - sphere.center.1, camera_origin.2 - sphere.center.2);



    todo!();
}
fn dot_product(a: &(i32,i32,i32), b: &(i32,i32,i32)) {
    let mut product = 0.0;
    for i in 0..2{
        product += a.i * b.i
    }
}

fn canvas_to_viewport(x: i32, y: i32, view_width: i32, view_height: i32, distance: i32) -> (i32, i32, i32) {

    let view_x = (x*view_width/(WIDTH as i32));
    let view_y = (y*view_height/(HEIGHT as i32));
    let view = (view_x, view_y, distance);
    return view;
    
}