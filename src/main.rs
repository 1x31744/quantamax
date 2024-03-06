extern crate sdl2;
// Add this import for the image crate
extern crate image;

use indicatif::{ProgressBar, ProgressStyle};

use std::sync::atomic::Ordering;
use rand::Rng;
use cond_utils::Between;
use sdl2::pixels;
use sdl2::pixels::PixelFormatEnum;
use sdl2::render;
use sdl2::render::Canvas;
use sdl2::render::TextureAccess;
use sdl2::render::TextureCreator;
use sdl2::sys::Atom;
use sdl2::video::WindowContext;
use std::clone;
use std::io::BufRead;
use std::num;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;
use std::vec;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use std::time::Duration;
use sdl2::rect::Rect;
use std::sync::mpsc::{self, Sender};
use sdl2::render::Texture;

const WIDTH: u32 = 1000;
const HEIGHT: u32 = 1000;
const UPDATE_FREQUENCY: u64 = 10000; // Update progress bar every 100 iterations
//TODO: implement a function that checks intersects on its own as this is done multiple times in the code
//TODO: make easy setting between certain features, like shadow_smoothing, global_illumination, right now shadow smoothing depends on global illumnation
#[derive(PartialEq, Clone)]
struct Sphere {
    pub radius: f64,
    pub center: [f64; 3],
    pub color: [u8; 3],
    pub specularity: f64,
    pub reflectivity: f64
}
impl Sphere {
    pub fn new(radius: f64, center: [f64; 3], color: [u8; 3], specularity: f64, reflectivity: f64) -> Sphere {
        Sphere {
            radius: radius,
            center: center,
            color: color,
            specularity: specularity,
            reflectivity: reflectivity
        }
    }
}

#[derive(PartialEq)]
struct Light {
    pub typ: String,
    pub intensity: f64,
    pub position: [f64; 3],
    pub direction: [f64; 3],
    pub radius: f64,
    pub sphere: Sphere,
}
impl Light {
    pub fn new (intensity: f64, typ: String, position: [f64; 3], direction: [f64; 3], radius: f64) -> Light {
        Light {
            intensity: intensity,
            typ: typ,
            position: position.clone(),
            direction: direction,
            radius: radius,
            sphere: Sphere::new(radius, position, [(225.0 * intensity) as u8,(225.0 * intensity) as u8,(225.0 * intensity) as u8], -1.0, 0.0)
        }
    }
}

pub fn main() {

    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("quantamax", WIDTH, HEIGHT)
        .position_centered()
        .build()
        .unwrap();

    //let res_multiplier: u32 = 1; //temp
    let mut resolution: [u32; 2] = [200, 200]; 
    let mut canvas = window.into_canvas().build().unwrap();
    //canvas.set_logical_size(resolution[0], resolution[1]).expect("could not set res");
    let mut event_pump = sdl_context.event_pump().unwrap();
    let texture_creator = canvas.texture_creator();
    let mut render_texture = texture_creator.create_texture_streaming(PixelFormatEnum::RGB24,resolution[0] + 1, resolution[1] + 1).map_err(|e| e.to_string()).unwrap();
    //let mut surface = window.surface(&event_pump);

    let window_width_half: i32 = (resolution[0]/2) as i32;
    let window_height_half: i32 = (resolution[1]/2) as i32;
    let sphere1 = Sphere::new(1.0, [0.0,-1.5, 1.5], [225,0,0] , -1.0, 0.5);
    let sphere2 = Sphere::new(1.0, [-1.5,-0.5,1.7], [255,192,203], 500.0, 0.5);
    let sphere3 = Sphere::new(1.0, [1.5,-0.5, 1.7], [0,0,225], 500.0, 0.5);
    let shere5 = Sphere::new(5000.0, [0.0,-5001.0,0.0], [225,225,0], -1.0, 0.2);
    let spheres = vec![sphere2,sphere1,sphere3,shere5];

    let light = Light::new(0.7, String::from("Point"), [-5.0,1.0,0.0], [0.0,0.0,0.0], 1.0);
    let light2 = Light::new(0.7, String::from("Point"), [5.0,1.0,0.0], [-1.0, -9.0, -2.0], 1.0); 
    let light3 = Light::new(0.2, String::from("Ambient"), [0.0,0.0,0.0], [0.0,0.0,0.0], 1.0);
    let lights: Vec<Light> = vec![light, light2, light3];

    let fov: f64 = 180 as f64;
    let view_width = (resolution[0] as f64) * fov;
    let view_height = (resolution[1] as f64) * fov;
    let mut camera_position = [0.0,0.0,-1.0];

    let mut camera_rotation_y: f64 = 0.0;
    let mut camera_rotation_x: f64 = 0.0;
    let camera_rotation_z: f64 = 0.0;

    let z_rotation_matrix: [[f64; 3]; 3] =[[f64::cos(camera_rotation_z).round(), -f64::sin(camera_rotation_z).round(), 0.0],
                                               [f64::sin(camera_rotation_z).round(), f64::cos(camera_rotation_z).round(), 0.0],
                                               [0.0,0.0,1.0]]; 
    let mut x_rotation_matrix: [[f64; 3]; 3] = [[1.0, 0.0, 0.0],
                                               [0.0, f64::cos(camera_rotation_x).round(), -f64::sin(camera_rotation_x).round()],
                                               [0.0, f64::sin(camera_rotation_x).round(), f64::cos(camera_rotation_x).round()]];
    let mut y_rotation_matrix: [[f64; 3]; 3] = [[f64::cos(camera_rotation_y).round(), 0.0, f64::sin(camera_rotation_y).round()],
                                                [0.0,1.0,0.0],
                                                [-f64::sin(camera_rotation_y).round(), 0.0, f64::cos(camera_rotation_y).round()]];


    //config that changes in runtime
    let mut r_pressed: bool = true;
    let mut render_chance = 200;
    let mut global_illumination = false;
    let mut smooth_shadows = false;
    let mut r_change_safety = true;
    'running: loop {
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..} |
                Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                    break 'running
                },
                Event::KeyDown {keycode: Some(Keycode::R), ..} => {
                    if r_change_safety {
                        if r_pressed {
                            render_chance = 200;
                            global_illumination = false;
                            smooth_shadows = false;
                            r_pressed = false;
                        }
                        else if r_pressed == false {
                            resolution = [800,800];
                            render_texture = texture_creator.create_texture_streaming(PixelFormatEnum::RGB24,resolution[0] + 1, resolution[1] + 1).map_err(|e| e.to_string()).unwrap();
                            render_chance = 1000;
                            global_illumination = true;
                            smooth_shadows = true;
                            r_pressed = true;
                        }
                    }
                    r_change_safety = false;
                },
                Event::KeyUp {keycode: Some(Keycode::R), ..} => {
                    r_change_safety = true;
                },
                Event::KeyDown {keycode: Some(Keycode::W), ..} => {
                    let original_move_direction = [0.0,0.0,0.1];
                    let with_x_rotation = matrix_multiplication(&x_rotation_matrix, original_move_direction);
                    let with_y_rotation = matrix_multiplication(&y_rotation_matrix, with_x_rotation);
                    camera_position = vec3_addition(&with_y_rotation, &camera_position);
                },
                Event::KeyDown {keycode: Some(Keycode::S), ..} => {
                    let original_move_direction = [0.0,0.0,-0.1];
                    let with_x_rotation = matrix_multiplication(&x_rotation_matrix, original_move_direction);
                    let with_y_rotation = matrix_multiplication(&y_rotation_matrix, with_x_rotation);
                    camera_position = vec3_addition(&with_y_rotation, &camera_position);
                },
                Event::KeyDown {keycode: Some(Keycode::D), ..} => {
                    let original_move_direction = [0.1,0.0,0.0];
                    let with_x_rotation = matrix_multiplication(&x_rotation_matrix, original_move_direction);
                    let with_y_rotation = matrix_multiplication(&y_rotation_matrix, with_x_rotation);
                    camera_position = vec3_addition(&with_y_rotation, &camera_position);
                },
                Event::KeyDown {keycode: Some(Keycode::A), ..} => {
                    let original_move_direction = [-0.1,0.0,0.0];
                    let with_x_rotation = matrix_multiplication(&x_rotation_matrix, original_move_direction);
                    let with_y_rotation = matrix_multiplication(&y_rotation_matrix, with_x_rotation);
                    camera_position = vec3_addition(&with_y_rotation, &camera_position);
                },
                Event::KeyDown {keycode: Some(Keycode::Space), ..} => {
                    let original_move_direction = [0.0,0.1,0.0];
                    let with_x_rotation = matrix_multiplication(&x_rotation_matrix, original_move_direction);
                    let with_y_rotation = matrix_multiplication(&y_rotation_matrix, with_x_rotation);
                    camera_position = vec3_addition(&with_y_rotation, &camera_position);
                },
                Event::KeyDown {keycode: Some(Keycode::Z), ..} => {
                    let original_move_direction = [0.0,-0.1,0.0];
                    let with_x_rotation = matrix_multiplication(&x_rotation_matrix, original_move_direction);
                    let with_y_rotation = matrix_multiplication(&y_rotation_matrix, with_x_rotation);
                    camera_position = vec3_addition(&with_y_rotation, &camera_position);
                },
                Event::KeyDown {keycode: Some(Keycode::Up), ..} => {
                    camera_rotation_x -= 0.1;
                    x_rotation_matrix = [[1.0, 0.0, 0.0],
                                               [0.0, f64::cos(camera_rotation_x), -f64::sin(camera_rotation_x)],
                                               [0.0, f64::sin(camera_rotation_x), f64::cos(camera_rotation_x)]];
                } 
                Event::KeyDown {keycode: Some(Keycode::Down), ..} => {
                    camera_rotation_x += 0.1;
                    x_rotation_matrix = [[1.0, 0.0, 0.0],
                                               [0.0, f64::cos(camera_rotation_x), -f64::sin(camera_rotation_x)],
                                               [0.0, f64::sin(camera_rotation_x), f64::cos(camera_rotation_x)]];
                } 
                Event::KeyDown {keycode: Some(Keycode::Left), ..} => {
                    camera_rotation_y -= 0.1;
                    y_rotation_matrix = [[f64::cos(camera_rotation_y), 0.0, f64::sin(camera_rotation_y)],
                                                [0.0,1.0,0.0],
                                                [-f64::sin(camera_rotation_y), 0.0, f64::cos(camera_rotation_y)]];
                } 
                Event::KeyDown {keycode: Some(Keycode::Right), ..} => {
                    camera_rotation_y += 0.1;
                    y_rotation_matrix = [[f64::cos(camera_rotation_y), 0.0, f64::sin(camera_rotation_y)],
                                                [0.0,1.0,0.0],
                                                [-f64::sin(camera_rotation_y), 0.0, f64::cos(camera_rotation_y)]];
                } 


                _ => {}
            }
        }

        canvas.clear();

        // ! a pixel buffer can be an array of u8's, in order R, G, B. first 3 are (0,0)
        

        let multithreading = true;
        if !multithreading {
            render_texture.with_lock(None, |buffer: &mut [u8], pitch: usize| {
                for x in -window_width_half..window_width_half {
                    for y in -window_height_half..window_height_half {
                        let view = canvas_to_viewport(x as f64, y as f64, view_width, view_height as f64, 1.0, &z_rotation_matrix, &x_rotation_matrix, &y_rotation_matrix);
                        let color = trace_ray(&camera_position, view, 0.0, f64::INFINITY, &lights, &spheres, 3.0, 2.0, global_illumination, smooth_shadows);
                        let canvas_coords = transfer_coords(x, y, window_height_half, window_height_half);      
                        let offset = canvas_coords.1 as usize * pitch as usize + canvas_coords.0 as usize * 3;
                        buffer[offset] = color[0];
                        buffer[offset + 1] = color[1];
                        buffer[offset + 2] = color[2];
                   }
                }
            }).expect("what?"); 
            canvas.copy(&render_texture, None, Rect::new(0, 0, WIDTH, HEIGHT)).expect("could not draw to texture");
        }
        if multithreading {
            //create channels for communication between threads
            let (tx, rx) = mpsc::channel();

            // Create a shared integer variable wrapped in Arc and Mutex
            let counter = Arc::new(AtomicUsize::new(0));
            //let pb = Arc::new(Mutex::new(ProgressBar::new((resolution[0] * resolution[1]) as u64)));
            //pb.lock().unwrap().set_style(ProgressStyle::default_bar()
            //.template("[{elapsed_precise}] {bar:40.red/pink} {percent}% {pos}/{len} {msg}")
            //.progress_chars("#>-"));

            //spawn multiple threads to update different sections of the texture
            let num_of_threads = 10;
            let handles: Vec<_> = (0..num_of_threads)
            .map(|i| {
                let counter_clone = Arc::clone(&counter);
                //let pb_clone = Arc::clone(&pb);

                let tx = tx.clone();
                let render_chance = render_chance.clone();
                let mut camera_position = camera_position.clone();
                let z_rotation_matrix = z_rotation_matrix.clone();
                let x_rotation_matrix = x_rotation_matrix.clone();
                let y_rotation_matrix = y_rotation_matrix.clone();
                let global_illumination = global_illumination.clone();
                let smooth_shadows = global_illumination.clone();
                let res = resolution.clone();
                let thread_num = num_of_threads.clone();
                thread::spawn(move || {
                    update_texture_region(i, tx, render_chance, &mut camera_position, &z_rotation_matrix,
                    &x_rotation_matrix, &y_rotation_matrix, global_illumination, smooth_shadows, res, counter_clone, thread_num);
                })
            }).collect();

            for handle in handles {
                for event in event_pump.poll_iter() {
                    match event {
                        Event::Quit {..} |
                        Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                            break 'running
                        },
                        _ => {}
                    }
                }
                handle.join().unwrap();
            }
            //recieve updates from threads and apply to texture
            for update in rx.try_iter() {
                render_texture.with_lock(None, |buffer: &mut [u8], pitch: usize| {
                    for (x, y, color) in update {
                        let offset = y as usize * pitch as usize + x as usize * 3;
                        buffer[offset] = color.0;
                        buffer[offset + 1] = color.1;
                        buffer[offset + 2] = color.2;
                    }
                }).expect("why error?");
            } 

            // Save the texture as an image and apply fxxa
            if global_illumination {
                // ! apply fxxa first
                let mut buffer: Vec<u8> = vec![0; ((resolution[0] * resolution[1] * 3)) as usize]; // Assuming RGB format
                let mut pixel_pitch: usize = 0;
                render_texture.with_lock(None, |data, pitch: usize| {
                    // Populate the buffer with pixel data
                    for y in 0..resolution[1] {
                        for x in 0..resolution[0] {
                            let offset = (y * resolution[0] + x) as usize * 3;
                            let r = data[y as usize * pitch as usize + x as usize * 3];
                            let g = data[y as usize * pitch as usize + x as usize * 3 + 1];
                            let b = data[y as usize * pitch as usize + x as usize * 3 + 2];
                            buffer[offset] = r;
                            buffer[offset + 1] = g;
                            buffer[offset + 2] = b;
                        }
                    }
                    pixel_pitch = pitch;
                }).unwrap();
                let mut modified_render_texture: Texture<'_>  = texture_creator.create_texture_streaming(PixelFormatEnum::RGB24,resolution[0] + 1, resolution[1] + 1).map_err(|e| e.to_string()).unwrap();

                render_texture.with_lock(None, |pixels, pitch| {
                    modified_render_texture.update(None, pixels, pitch).expect("could not copy section");
                }).expect("could not copy section");

                buffer = apply_fxaa(&texture_creator, &mut canvas, &mut modified_render_texture, &mut buffer, &pixel_pitch);

                // ! make image
                let img = image::RgbImage::from_raw(resolution[0], resolution[1], buffer).unwrap();
                let path = "output.png";
                img.save(path).unwrap();
            }

            canvas.copy(&render_texture, None, Rect::new(0, 0, WIDTH, HEIGHT)).expect("could not draw");  
            canvas.present();

            for _ in 0..num_of_threads{
                tx.send(Vec::new()).unwrap();
            }
        }
        if global_illumination {
            println!("render finished")
        }
        if global_illumination == true {println!("render finshed, press esc to quit or r to go back into low render mode")}
        while global_illumination == true {
            for event in event_pump.poll_iter() {
                match event {
                    Event::Quit {..} |
                    Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                        break 'running
                    },
                    Event::KeyDown {keycode: Some(Keycode::R), ..} => {
                                render_chance = 200;
                                global_illumination = false;
                                smooth_shadows = false;
                                r_pressed = false;

                        break;
                    },
                    _ => {}
                }
            }
        }
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60));
        //break;
    }
}

fn apply_fxaa<'a>(
    texture_creator: &'a TextureCreator<WindowContext>,
    canvas: &mut Canvas<sdl2::video::Window>,
    texture: &mut Texture<'a>,
    buffer: &mut Vec<u8>,
    pitch: &usize
) -> Vec<u8>{
    // Extract the pixel data from the texture
    let query = texture.query();
    let (width, height) = (query.width as usize, query.height as usize);
    let mut pixels: Vec<u8> = vec![0; ((width * height * 3)) as usize];

    // Apply FXAA
    return fxaa(buffer, width, height, pitch);

    // Update the texture with the modified pixels

    //texture.update(None, &pixels, (width * 3) as usize * 3 as usize).map_err(|e| e.to_string()).expect("hmm these are some pretty good error");
    //canvas.copy(texture, None, Rect::new(0, 0, WIDTH, HEIGHT)).expect("could not draw");  

}
fn fxaa(pixels: &Vec<u8>, width: usize, height: usize, pitch: &usize) -> Vec<u8> {
    let mut luma_current;
    let mut luma_down: f32;
    let mut luma_right;
    let mut new_buffer = pixels.clone();
    
    let r_lum_mult = 0.2126;
    let g_lum_mult = 0.7152;
    let b_lum_mult = 0.0722;

    for y in 0..(height-3) {
        for x in 0..(width-6) {
            let offset = (y * width + x) as usize * 3; // RGB24 format
            let down_offset = ((y + 1) * width + x) as usize * 3;

            // Calculate luminance for the current pixel and its neighbors
            luma_current = r_lum_mult * pixels[offset] as f32 + g_lum_mult * pixels[offset + 1] as f32 + b_lum_mult * pixels[offset+2] as f32;
            if x > 0 {
                luma_right = r_lum_mult * pixels[offset + 3] as f32 + g_lum_mult * pixels[offset + 4] as f32 + b_lum_mult * pixels[offset + 5] as f32;
            } else {
                luma_right = luma_current;
            }
            if y > 0 {
                luma_down = r_lum_mult * pixels[down_offset] as f32 + g_lum_mult * pixels[down_offset + 1] as f32 + b_lum_mult * pixels[down_offset + 2] as f32;
            } else {
                luma_down = luma_current;
            }

            //println!("Current pixel: ({}, {})", x, y);
            //println!("Luminance current: {}", luma_current);
            //println!("Luminance down: {}", luma_down);
            //println!("Luminance right: {}", luma_right);
    

            // Calculate contrast between the current pixel and its neighbors
            //println!("{}, {}", luma_current, luma_right);
            let contrast_down = (luma_current - luma_down).abs();
            let contrast_right = (luma_current - luma_right).abs();
            //println!("contrast down is: {}", contrast_down);
            //println!("contrast right is: {}", contrast_right);
            //println!("{}",contrast_right);
            
            // Apply FXAA if the contrast is high
            if contrast_right >= 120.0 || (luma_right == 0.0 && luma_current > 30.0) || (luma_current == 0.0 && luma_right > 30.0){ // TODO: try doing this seperately, handle right contrast differently to down, otherwise bad things considered
                // Apply a simple blur filter to the pixel;
                for i in 0..3 {
                    let neighbor_pixel_right = pixels[offset + 3 + i];
                    new_buffer[offset + i] = (pixels[offset + i] + neighbor_pixel_right) / 2;
                }
                println!("this happen")
            }
            if contrast_down >= 120.0 || (luma_down == 0.0 && luma_current > 30.0) || (luma_current == 0.0 && luma_down > 30.0) {
                for i in 0..3 {
                    let neighbor_pixel_down = pixels[down_offset + i];
                    new_buffer[offset + i] = (pixels[offset + i] + neighbor_pixel_down) / 2;
                }
                println!("this happen")

            // ! i think the issue might be that changes in the average color are effecting right and below calculations, therefore the we change should be different
            // ! from the one we analyse
            }
        }
    }
    return new_buffer.to_vec()
}

fn update_texture_region(thread_id: usize, tx: Sender<Vec<(u32, u32, (u8, u8, u8))>>, render_chance_max: i32, camera_position: &mut [f64; 3], z_rotation_matrix: &[[f64; 3]; 3], x_rotation_matrix: &[[f64; 3]; 3]
, y_rotation_matrix: &[[f64; 3]; 3], global_illumination: bool, smooth_shadows: bool, resolution: [u32; 2], percent_counter: Arc<AtomicUsize>, num_of_threads: usize){
    // --temp--
    //define variables

    let res_width_half: i32 = (resolution[0]/2) as i32;
    let res_height_half: i32 = (resolution[1]/2) as i32;
    let sphere1 = Sphere::new(1.0, [0.0,-1.5, 1.5], [225,0,0] , -1.0, 0.5);
    let sphere2 = Sphere::new(1.0, [-1.5,-0.5,1.7], [255,192,203], 500.0, 0.5);
    let sphere3 = Sphere::new(1.0, [1.5,-0.5, 1.7], [0,0,225], 500.0, 0.5);
    let shere5 = Sphere::new(5000.0, [0.0,-5001.0,0.0], [225,225,0], -1.0, 0.2);
    let spheres = vec![sphere2,sphere1,sphere3,shere5];

    let light = Light::new(0.7, String::from("Point"), [-5.0,1.0,0.0], [0.0,0.0,0.0], 0.5);
    let light2 = Light::new(0.7, String::from("Point"), [5.0,1.0,0.0], [-1.0, -9.0, -2.0], 0.5); 
    let lights: Vec<Light> = vec![light, light2];

    let fov: f64 = 180 as f64;
    let view_width = (resolution[0] as f64) * fov;
    let view_height = (resolution[1] as f64) * fov;

    let num_of_threads = num_of_threads as u32;
    let render_width = resolution[0]/num_of_threads; //width/numofthreads


    let mut update: Vec<(u32, u32, (u8, u8, u8))> = Vec::new();
    let starting_point = -res_width_half + (render_width * (thread_id) as u32) as i32;

    let max_negate_id = (num_of_threads - (thread_id + 1) as u32) as i32;
    let end_point = res_width_half - (render_width * (max_negate_id) as u32) as i32;

    
    let max_pixels = resolution[0] * resolution[1];
    let mut rng = rand::thread_rng();
    for x in starting_point..end_point {
        for y in -res_height_half..res_height_half {
            let color: [u8; 3];
            let render_chance = rng.gen_range(0..1000);
            if render_chance < render_chance_max {
                let view = canvas_to_viewport(x as f64, y as f64, view_width, view_height as f64, 1.0, &z_rotation_matrix, &x_rotation_matrix, &y_rotation_matrix);
                color = trace_ray(&camera_position, view, 0.0, f64::INFINITY, &lights, &spheres, 3.0, 1.0, global_illumination, smooth_shadows);
            }
            else  {
                color = [0,0,0];
            }
            let canvas_coords = transfer_coords(x, y, res_width_half, res_height_half);
            percent_counter.fetch_add(1, Ordering::SeqCst); // add 1
            let value = percent_counter.load(Ordering::SeqCst);
            if (value as u64 % UPDATE_FREQUENCY == 0 || value == 1) && global_illumination{
                println!("{}%", (value as f32/max_pixels as f32) * 100.0);
            }

            update.push((canvas_coords.0 as u32, canvas_coords.1 as u32, (color[0], color[1], color[2])));
       }
    }
    tx.send(update).unwrap();

        
}

fn transfer_coords(x: i32, y: i32, res_width_half: i32, res_height_half: i32) -> (i32, i32) {
    let canvas_x: i32 = (x + res_width_half) as i32;
    let canvas_y: i32 = (res_height_half - y) as i32;
    return (canvas_x, canvas_y)
}

fn trace_ray(camera_origin: &[f64; 3], view: [f64; 3], tmin: f64, tmax: f64, lights: &Vec<Light>, spheres: &Vec<Sphere>, recursion_depth_reflection: f64, recursion_depth_indirect: f64
, global_illumination: bool, smooth_shadows: bool) -> [u8; 3] {
    let mut closest_sphere: Option<&Sphere> = None;
    let mut closest_t: f64 = f64::INFINITY;
    let mut intersects: [f64; 2];
    let mut is_light_source: bool = false;

    for sphere in spheres{
        intersects = intersect_ray_sphere(&camera_origin, &view, sphere); //two intersections
        if intersects[0].within(tmin, tmax) && intersects[0] < closest_t {
            closest_t = intersects[0];
            closest_sphere = Some(sphere);
            is_light_source = false;
        }
        if intersects[1].within(tmin, tmax) && intersects[1] < closest_t {
            closest_t = intersects[1];
            closest_sphere = Some(sphere);
            is_light_source = false;
        }
    }
    
    //do the same but for light sources
    for light in lights {
        intersects = intersect_ray_sphere(&camera_origin, &view, &light.sphere);
        if intersects[0].within(tmin, tmax) && intersects[0] < closest_t {
            closest_t = intersects[0];
            closest_sphere = Some(&light.sphere);
            is_light_source = true;
        }
        if intersects[1].within(tmin, tmax) && intersects[1] < closest_t {
            closest_t = intersects[1];
            closest_sphere = Some(&light.sphere);
            is_light_source = true;
        }
    }

    if closest_sphere == None {
        return [0, 0, 0]
    }
    //assert_eq!(intersects,vec![0.0, 0.0]);

    let unwraped_sphere = closest_sphere.unwrap();

    if is_light_source {
        return unwraped_sphere.color.clone();
    }

    //--values for light intensity calculation--
    let intersection: [f64; 3] = vec3_addition(&camera_origin, &vec3_multiply_by_float(&view, closest_t));
    let mut normal = vec3_negation(&intersection, &unwraped_sphere.center);
    normal = vec3_divide_by_float(&normal, vec3_length(&normal));
    //--values for light intensity calculation--

    let light_intensity = compute_lighting(&intersection, &normal, lights, vec3_multiply_by_float(&view, -1.0), unwraped_sphere.specularity, spheres, smooth_shadows);
    
    let local_color = multiply_color_by_float(&unwraped_sphere.color, light_intensity);

    //compute reflection
    let reflectivity = unwraped_sphere.reflectivity;

    if recursion_depth_reflection < 0.0 || reflectivity == 0.0 {
        return local_color
    }

    //pythagoras
    let opposite_ray = vec3_multiply_by_float(&view, -1.0);
    let reflection_with_view = vec3_multiply_by_float(&normal, 2.0 * dot_product(&normal, &opposite_ray));
    let reflect_ray = vec3_negation(&reflection_with_view, &opposite_ray);


    let reflected_color = trace_ray(&intersection, reflect_ray, 0.001, f64::INFINITY, lights, spheres, recursion_depth_reflection-1.0, recursion_depth_indirect, global_illumination, smooth_shadows);
    
    let mut indirect_lighting = [0,0,0];
    if global_illumination == true {
        //indirect lighting comes after reflectivity? (also recursive)
        let num_of_samples = 40.0;
        let mut rng = rand::thread_rng();
        if recursion_depth_indirect > 0.0 { //use indirect's own recursion depth
            for _ in 0..num_of_samples as u32 { //num of indirect samples
                //produce random direction using the unit circle
                let random_direction = [rng.gen::<f64>() * 2.0 - 1.0, rng.gen::<f64>() * 2.0 - 1.0, rng.gen::<f64>() * 2.0 - 1.0];
                let normalized_random_direction = normalize(&random_direction);
                let indirect_color = trace_ray(&intersection, normalized_random_direction , 0.001, f64::INFINITY, lights, spheres, 0.0, recursion_depth_indirect - 1.0, global_illumination, smooth_shadows);
                for i in 0..indirect_color.len() {
                    indirect_lighting[i] = indirect_color[i] + indirect_lighting[i] //add indirect color to the total indirect ligting
                }
            }
            //divide indirect by the number of samples
            indirect_lighting = multiply_color_by_float(&indirect_lighting, 1.0/num_of_samples);
        }
        else if recursion_depth_indirect < 0.0 {return indirect_lighting;}
    }
    //factor in reflectivity and indirect lighting
    let real_local = multiply_color_by_float(&local_color, 1.0-reflectivity);
    let real_reflective = multiply_color_by_float(&reflected_color, reflectivity);
    let mut global_color = [0,0,0];
    for i in 0..real_local.len() {
        global_color[i] = real_local[i] + real_reflective[i];
        if global_illumination == true {
            global_color[i] = global_color[i] + indirect_lighting[i];
        }
        //println!("{:?}", global_color);
        if global_color[i] > 225 { global_color[i] = 225}
    }
    return global_color;
}

// ! may use later, as finding the intersects is used in alot of places in code, using a function for it may be a little better, especially when adding other mediums
/*  
fn find_intersects(camera_origin: &Vec<f64>, view: &Vec<f64>, tmin: f64, tmax: f64, spheres: & mut Vec<Sphere>) -> (Option<&mut Sphere>, f64) {
    let mut closest_sphere: Option<& mut Sphere> = None;
    let mut closest_t: f64 = f64::INFINITY;
    let mut intersects: Vec<f64>;
    for sphere in spheres{
        intersects = intersect_ray_sphere(&camera_origin, &view, sphere); //two intersections
        if intersects[1].within(tmin, tmax) && intersects[1] < closest_t {
            closest_t = intersects[1];
            closest_sphere = Some(sphere)
        }
        if intersects[0].within(tmin, tmax) && intersects[0] < closest_t {
            closest_t = intersects[0];
        }
    }
    return  (closest_sphere, closest_t);
}
*/

fn intersect_ray_sphere(camera_origin: &[f64; 3], view: &[f64; 3], sphere: &Sphere) -> [f64; 2] {
    let radius = sphere.radius;
    let camera_to_center = [camera_origin[0] - sphere.center[0], camera_origin[1] - sphere.center[1], camera_origin[2] - sphere.center[2]];
    //let camera_to_center = vec![sphere.center[0] - camera_origin[0], sphere.center[1] - camera_origin[1], sphere.center[2] - camera_origin[2]];
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
    let a: f64 = dot_product(&view, &view);
    let b :f64 = 2.0 * dot_product(&view, &camera_to_center);
    let c: f64 = (dot_product(&camera_to_center, &camera_to_center)) - radius.powi(2);

    let discriminant: f64 = b.powi(2)  - (4.0*a*c);
    if discriminant < 0.0 {
        return [f64::INFINITY, f64::INFINITY];
    }
    //then do quadratic formula
    let t1: f64 = (-b + f64::sqrt(discriminant)) / (2.0*a);
    let t2: f64 = (-b - f64::sqrt(discriminant)) / (2.0*a);
    return [t1, t2];
}


fn canvas_to_viewport(x: f64, y: f64, view_width: f64, view_height: f64, distance: f64, z_rotation_matrix: &[[f64; 3]; 3], x_rotation_matrix: &[[f64; 3]; 3], y_rotation_matrix: &[[f64; 3]; 3]) -> [f64; 3] {

    let view_x:f64 = (x*(WIDTH as f64)/view_width) as f64;
    let view_y:f64 = (y*(HEIGHT as f64)/view_height) as f64;
    let mut view = [view_x, view_y, distance];

    view = matrix_multiplication(z_rotation_matrix, view);
    view = matrix_multiplication(x_rotation_matrix, view);
    view = matrix_multiplication(y_rotation_matrix, view);

    return view;
    
}

fn compute_lighting(intersection: &[f64; 3], normal: &[f64; 3], light: &Vec<Light>, to_cam: [f64; 3], specularity: f64, spheres: &Vec<Sphere>, shadow_smoothing: bool) -> f64 {
    let mut i: f64 = 0.0;
    let mut light_direction: [f64; 3];
    let mut t_max: f64;
    for light in light {
        if light.typ == "Ambient" {
            i += light.intensity;
        }
        else {
            if light.typ == "Point" {
                light_direction = [light.position[0] - intersection[0], light.position[1] - intersection[1], light.position[2] - intersection[2]];
                t_max = 1.0;
            }
            else {
                light_direction = light.direction.clone();
                t_max = f64::INFINITY
            }
            if !shadow_smoothing { // ! old implementation for performance
                let mut shadow_sphere: Option<&Sphere> = None;
                let mut shadow_t: f64 = f64::INFINITY;
                let mut intersects;
                //TODO: there would be multiple light directions;
                for sphere in spheres{
                    intersects = intersect_ray_sphere(&intersection, &light_direction, sphere); //two intersections
                    if intersects[0].within(0.001, t_max) && intersects[0] < shadow_t {
                        shadow_t = intersects[0];
                        shadow_sphere = Some(sphere);
                    }
                    if intersects[1].within(0.001, t_max) && intersects[1] < shadow_t {
                        shadow_t = intersects[1];
                        shadow_sphere = Some(sphere)
                    }
                }
                if !(shadow_sphere == None) {
                    let normal_dot_light = dot_product(&normal, &light_direction);
                    if normal_dot_light > 0.0 {
                        //this is true because of trigonometry
                        //cos = distance to cam / normal + light direction
                        i += (light.intensity * normal_dot_light)/((vec3_length(&normal) * vec3_length(&light_direction)))
                    }
                    //for specular reflection
                    if specularity != -1.0 {
                        let reflection =  vec3_negation(&vec3_multiply_by_float(&normal, dot_product(&normal, &light_direction) * 2.0), &light_direction);
                        let reflection_cam_distance = dot_product(&reflection, &to_cam);
                        if reflection_cam_distance > 0.0 {
                            i += light.intensity * f64::powf(reflection_cam_distance/(vec3_length(&reflection)* vec3_length(&to_cam)), specularity)
                        }
                    }
                }
            }
            else { // ! shadow smoothing implementation
                let light_percent = 1.0 - compute_shadows_percentage(light, spheres, &light.position,  intersection);    
                //for diffuse reflection
                let normal_dot_light = dot_product(&normal, &light_direction);
                if normal_dot_light > 0.0 {
                    //this is true because of trigonometry
                    //cos = distance to cam / normal + light direction
                    i += ((light.intensity * normal_dot_light)/((vec3_length(&normal) * vec3_length(&light_direction)))) * light_percent as f64
                }
                //for specular reflection
                if specularity != -1.0 {
                    let reflection =  vec3_negation(&vec3_multiply_by_float(&normal, dot_product(&normal, &light_direction) * 2.0), &light_direction);
                    let reflection_cam_distance = dot_product(&reflection, &to_cam);
                    if reflection_cam_distance > 0.0 {
                        i += light.intensity * f64::powf(reflection_cam_distance/(vec3_length(&reflection)* vec3_length(&to_cam)), specularity)
                    }
                }
            }
        }
    }
    return i;
}
fn compute_shadows_percentage(light: &Light, spheres: &Vec<Sphere>, light_position: &[f64; 3], intersection: &[f64; 3])  -> f64{ //u8 is percentage illumination
    // ! compute angle from point, light_center and light edge
    let num_of_samples = 100;
    //let mut rng = rand::thread_rng();
    let mut successful_samples = 0;
    for _ in 0..num_of_samples {

        //sample point would be a random point on the circle
        let sample_point_x = light_position[0] + ((rand::random::<f64>() - 0.5) * 2.0) * light.radius;
        let sample_point_y = light_position[1] + ((rand::random::<f64>() - 0.5) * 2.0) * light.radius;
        let sample_point_z = light_position[2] + ((rand::random::<f64>() - 0.5) * 2.0) * light.radius;
        let sample_point = [sample_point_x, sample_point_y, sample_point_z];

        let shadow_ray = vec3_negation(&sample_point, &intersection);

        //set up variables for calculating intersections
        let mut shadow_sphere: Option<&Sphere> = None;
        let mut shadow_t: f64 = f64::INFINITY;
        let mut intersects;
        let t_max = 1.0;

        for sphere in spheres{
            intersects = intersect_ray_sphere(&intersection, &shadow_ray, sphere); //two intersections
            if intersects[0].within(0.001, t_max) && intersects[0] < shadow_t {
                shadow_t = intersects[0];
                shadow_sphere = Some(sphere);
            }
            if intersects[1].within(0.001, t_max) && intersects[1] < shadow_t {
                shadow_t = intersects[1];
                shadow_sphere = Some(sphere)
            }
        }
        if !(shadow_sphere == None)  {
            successful_samples += 1;
        }
    }
    return successful_samples as f64/num_of_samples as f64; 
} 
fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let mut product: f64 = 0.0;
    for i in 0..3{
        product = product +  (a[i] * b[i]);
    }
    return product;
}
fn normalize(a: &[f64; 3]) -> [f64; 3] {
    let normalized_scale = 1.0/(f64::sqrt(a[0].powi(2) + a[1].powi(2) + a[2].powi(2)));
    let normalized = vec3_multiply_by_float(&a, normalized_scale);
    return normalized;
}
//vector manipulation functions
fn vec3_addition(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    let mut product: [f64; 3] = [0.0,0.0,0.0];
    for i in 0..a.len() {
        product[i] = a[i] + b[i]
    }
    return product;
}
fn vec3_negation(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    let mut product: [f64; 3] = [0.0,0.0,0.0];
    for i in 0..a.len() {
        product[i] = a[i] - b[i]
    }
    return product;
}
fn vec3_multiply_by_float(vec3: &[f64; 3], multiplier: f64) -> [f64; 3] {
    let mut product: [f64; 3] = [0.0,0.0,0.0];
    for i in 0..vec3.len() {
        product[i] = multiplier * vec3[i]
    }
    return product
}
fn vec3_divide_by_float(vec3: &[f64; 3], multiplier: f64) -> [f64; 3] {
    let mut product: [f64; 3] = [0.0,0.0,0.0];
    for i in 0..vec3.len() {
        product[i] = vec3[i] / multiplier
    }
    return product

}
fn vec3_length(vec3: &[f64; 3]) -> f64 {
    return f64::sqrt(vec3[0].powi(2) + vec3[1].powi(2)+ vec3[2].powi(2))
}
fn multiply_color_by_float(vec3: &[u8; 3], multiplier: f64) -> [u8; 3] {
    let mut product: [u8; 3] = [0,0,0];
    for i in 0..vec3.len() {
        product[i] = ((multiplier) * vec3[i] as f64).round() as u8
    }
    return product
}

fn matrix_multiplication(transformation: &[[f64; 3]; 3], matrix: [f64; 3]) -> [f64; 3] {
    let mut new = [0.0,0.0,0.0];
    for i in 0..transformation.len() {
        let trans = &transformation[i];
        let product = dot_product(trans, &matrix);
        new[i] = product
    }
    return new
}