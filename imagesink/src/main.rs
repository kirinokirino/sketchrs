use glam::UVec2;
use speedy2d::color::Color;
use speedy2d::image::{ImageDataType, ImageHandle, ImageSmoothingMode};
use speedy2d::window::{WindowCreationOptions, WindowHandler, WindowHelper, WindowSize};
use speedy2d::{Graphics2D, Window};

use memmap2::Mmap;
use std::env;
use std::fs::File;
use std::io::Write;

const WIDTH: u16 = 640;
const HEIGHT: u16 = 480;

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let (width, height) = if args.len() == 2 {
        (
            args[0].parse::<u16>().unwrap_or(WIDTH),
            args[1].parse::<u16>().unwrap_or(HEIGHT),
        )
    } else {
        (WIDTH, HEIGHT)
    };
    let mut file = File::options()
        .create(true)
        .truncate(true)
        .read(true)
        .write(true)
        .open("/tmp/imagesink")
        .unwrap();
    let _ = file.write_all(vec![0; usize::from(width) * usize::from(height) * 4].as_slice());

    let mmap = unsafe { Mmap::map(&file).unwrap() };
    let _ = mmap.lock();

    let size = UVec2::new(u32::from(width), u32::from(height));
    let window_size = WindowSize::PhysicalPixels(size);
    let window_options = WindowCreationOptions::new_windowed(window_size, None)
        .with_decorations(false)
        .with_transparent(true)
        .with_resizable(false);
    let window = Window::new_with_options("FLOATING", window_options).unwrap();
    window.run_loop(MyWindowHandler {
        mmap,
        width,
        height,
    });
}

struct MyWindowHandler {
    //image: Option<ImageHandle>,
    mmap: Mmap,
    width: u16,
    height: u16,
}

impl WindowHandler for MyWindowHandler {
    fn on_draw(&mut self, helper: &mut WindowHelper, graphics: &mut Graphics2D) {
        let image: ImageHandle = graphics
            .create_image_from_raw_pixels(
                ImageDataType::RGBA,
                ImageSmoothingMode::NearestNeighbor,
                UVec2::new(u32::from(self.width), u32::from(self.height)),
                &self.mmap[..],
            )
            .unwrap();

        helper.set_size_pixels(*image.size());

        graphics.clear_screen(Color::TRANSPARENT);

        graphics.draw_image((0.0, 0.0), &image);
        helper.request_redraw();
    }

    // If desired, on_mouse_move(), on_key_down(), etc...
}
