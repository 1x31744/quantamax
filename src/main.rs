extern crate sdl2;

use cond_utils::Between;
use sdl2::pixels::Color;
use std::fs::File;
use std::io::prelude::*;
use std::mem::Discriminant;
use std::vec;
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
    pub radius: i32,
    pub center: Vec<i32>,
    pub color: Vec<u8>,

}
impl Sphere {
    pub fn new(radius: i32, center: Vec<i32>, color: Vec<u8>) -> Sphere {
        Sphere {
            radius: radius,
            center: center,
            color: color,
        }
    }
}


fn put_pixel(canvas: &mut Canvas<Window>,x: i32, y:i32, color: Vec<u8>) {
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
    canvas.set_draw_color(Color::RGB(color[0], color[1], color[2]));
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

    let window_width_half: i32 = (WIDTH/2) as i32;
    let window_height_half: i32 = (HEIGHT/2) as i32;
    //let mut sphere1 = Sphere::new(100, vec![2,0, -5], vec![225,0,0]);
    let mut sphere2 = Sphere::new(2, vec![0,0,5], vec![0,225,0]);
    let mut spheres = vec![sphere2];
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

        //canvas.set_draw_color(Color::RGB(225, 225, 225));
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
                let view = canvas_to_viewport(x, y, (WIDTH*5) as i32, (HEIGHT*5) as i32, 1);
                //let debug_text = format!("x, y = {}, {} \n view = {:?} \n\n", x, y, viewport);
                //write!(debug, "{}", debug_text);
                let color = trace_ray(vec![0,0,0], view, 2, 200000000, &mut spheres);
                put_pixel(&mut canvas, x, y, color)
           }
           //write!(debug, "{}", "new line");
        }
        println!("done");
        //put_pixel(&mut canvas, 1, 1, vec![225, 0,0]);
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60));
        break;
    }
    ::std::thread::sleep(Duration::new(10, 1_000_000_000u32 / 60));
}

fn trace_ray(camera_origin: Vec<i32>, view: Vec<i32>, tmin: i32, tmax: i32, spheres: &mut Vec<Sphere>) -> Vec<u8> {
    let mut closest_sphere: Option<&mut Sphere> = None;
    let mut closest_t: i32 = 2147483647;
    for sphere in spheres{
        let intersects : Vec<f32> = intersect_ray_sphere(&camera_origin, &view, sphere); //two intersections
        if intersects[0].within(tmin as f32, tmax as f32) && intersects[0] < closest_t as f32 {
            closest_t = intersects[0] as i32;
            closest_sphere = Some(sphere);
        }
        if intersects[1].within(tmin as f32, tmax as f32) && intersects[1] < closest_t as f32 {
            closest_t = intersects[1] as i32;
            closest_sphere = Some(sphere)
        }
        if closest_sphere == None {
            return vec![225, 225, 225]
        }
    }
    //assert_eq!(intersects,vec![0.0, 0.0]);
    return closest_sphere.unwrap().color.clone();
}

fn intersect_ray_sphere(camera_origin: &Vec<i32>, view: &Vec<i32>, sphere: &mut Sphere) -> Vec<f32> {
    let radius = sphere.radius;
    //let camera_to_center = vec![camera_origin[0] - sphere.center[0], camera_origin[1] - sphere.center[1], camera_origin[2] - sphere.center[2]];
    let camera_to_center = vec![sphere.center[0] - camera_origin[0], sphere.center[1] - camera_origin[1], sphere.center[2] - camera_origin[2]];
    //println!("{:?}",camera_to_center);
    //ray equation moddeled by point = t*view + camera_origin
    //circle equation moddled by dot_product(point - center, point - center) = r^2
    //sub ray into circle: dot_product((t*view+camera_origin) - center, (t*view+camera_origin) - center) = r^2
    // cam_to center = camorgigin - center
    //plug in cam to center dot_product(camera_to_center + t*view, camera +t*view) = r^2
    //expand dot product: dot_product((cam_to_center + t*view), view) + dot_product((cam_to_center + t*view), t*view)
    //expand again: dot(cam_to_center, cam_to_center) + dot(t*view, cam_to_center) + dot(cam_to_center, t*view) + dot(t*view, t*view)
    //finally: t^2*dot(view,view) + t*2dot(view, am_to_center) + dot(cam_to_center, cam_to_center) - r^2 = 0
    //DO QUADRATIC, SOLVE FOR T
    //println!("{:?}, {:?}", camera_to_center, view);
    let a: f32 = dot_product(&view, &view) as f32;
    let b :f32 = 2.0 * dot_product(&view, &camera_to_center) as f32;
    let c: f32 = (dot_product(&camera_to_center, &camera_to_center) - radius.pow(2)) as f32;

    let discriminant: f32 = b.powf(2.0)  - 4.0*a*c;
    if discriminant < 0.0 {
        return vec![34359738368.0, 34359738368.0];
    }
    println!("{:?}", discriminant);
    //then do quadratic formula
    let t1: f32 = (-b + f32::sqrt(discriminant)) / (2.0*a);
    let t2: f32 = (-b - f32::sqrt(discriminant)) / (2.0*a);
    //print!("{:?}", t1);
    //println!("{}", t1);
    return vec![t1, t2];
}
fn dot_product(a: &Vec<i32>, b: &Vec<i32>) -> i32 {
    let mut product: i32 = 0;
    for i in 0..a.len(){
        product += a[i] * b[i]
    }
    return product;
}

fn canvas_to_viewport(x: i32, y: i32, view_width: i32, view_height: i32, distance: i32) -> Vec<i32> {

    let view_x = (x*view_width/(WIDTH as i32));
    let view_y = (y*view_height/(HEIGHT as i32));
    let view = vec![view_x, view_y, distance];
    //println!("{:?}", view);
    return view;
    
}