extern crate sdl2;

use sdl2::pixels::Color;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::render::Canvas;
use sdl2::video::Window;
use core::{f32, panic};
use std::time::Duration;
use sdl2::rect::Point;

const WIDTH: u32 = 1000;
const HEIGHT: u32 = 500;

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

pub fn main() {
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("quantamax", WIDTH, HEIGHT)
        .position_centered()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();
    let mut event_pump = sdl_context.event_pump().unwrap();
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
        println!("done");
        put_pixel(&mut canvas, 1, 1, (225,0,0));
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60));
    }
}