use glam::UVec2;
use memmap2::Mmap;
use std::env;
use std::fs::File;
use std::io::Write;
use std::time::{Duration, Instant};

use sdl2::event::Event;
use sdl2::pixels::PixelFormatEnum;
use sdl2::rect::Rect;

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

    let mmap_size = usize::from(width) * usize::from(height) * 4;
    let mut file = File::options()
        .create(true)
        .truncate(true)
        .read(true)
        .write(true)
        .open("/tmp/imagesink")
        .unwrap();
    let _ = file.write_all(&vec![0; mmap_size]);

    let mmap = unsafe { Mmap::map(&file).unwrap() };

    let sdl = sdl2::init().unwrap();
    let video = sdl.video().unwrap();
    let window = video
        .window("imagesink", width.into(), height.into())
        .position_centered()
        .opengl()
        .resizable()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().accelerated().present_vsync().build().unwrap();
    let texture_creator = canvas.texture_creator();
    let mut texture = texture_creator
        .create_texture_streaming(PixelFormatEnum::RGBA8888, width.into(), height.into())
        .unwrap();

    let mut event_pump = sdl.event_pump().unwrap();
    let mut last_frame = Instant::now();

    'running: loop {
        for event in event_pump.poll_iter() {
            if let Event::Quit { .. } = event {
                break 'running;
            }
        }

        // Lock and update texture from mmap
        texture
            .with_lock(None, |buffer: &mut [u8], _pitch: usize| {
                buffer.copy_from_slice(&mmap[..]);
            })
            .unwrap();

        canvas.clear();
        canvas.copy(&texture, None, Some(Rect::new(0, 0, width.into(), height.into()))).unwrap();
        canvas.present();

        // Frame cap (optional)
        let frame_time = last_frame.elapsed();
        if frame_time < Duration::from_millis(16) {
            std::thread::sleep(Duration::from_millis(16) - frame_time);
        }
        last_frame = Instant::now();
    }
}
