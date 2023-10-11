use std::io::Write;
use std::process::{Command, Stdio};

const RECORD: bool = false;
const WIDTH: usize = 640;
const HEIGHT: usize = 480;

use glam::Vec2;

fn main() {
    let ffmpeg_command = "/usr/bin/ffmpeg";
    let args = "-y -f rawvideo -vcodec rawvideo -s 640x480 -pix_fmt rgba -r 30 -i - -an -vcodec h264 -pix_fmt yuv420p -crf 15 /home/kirinokirino/Media/video.mp4".split(' ');
    let mut ffmpeg = if RECORD {
        Some(
            Command::new(ffmpeg_command)
                .args(args)
                .stdin(Stdio::piped())
                .spawn()
                .expect("failed to execute process")
                .stdin
                .take()
                .unwrap(),
        )
    } else {
        None
    };

    let mut palette: Vec<_> = ["#df8c00", "#d7ac64", "#36241e", "#4a3d35", "#747769"]
        .iter()
        .map(hex_to_rgb)
        .collect();
    palette.extend([[0, 0, 0, 0]].repeat(30));
    let mut canvas = Canvas::new(palette);
    //canvas.random();
    loop {
        canvas.update();
        canvas.display();
        if RECORD {
            ffmpeg.as_mut().map(|mut ffmpeg| {
                ffmpeg.write_all(
                    &canvas
                        .buffer.as_slice()
                        // .iter()
                        // .enumerate()
                        // .filter(|(i, v)| i % 4 != 3)
                        // .map(|(i, v)| *v)
                        // .collect::<Vec<u8>>(),
                )
            });
        }
        std::thread::sleep(std::time::Duration::from_millis(1));
    }
}

struct Canvas {
    pub buffer: Vec<u8>,
    palette: Vec<[u8; 4]>,
    pen_color: [u8; 4],
}

impl Canvas {
    pub fn new(palette: Vec<[u8; 4]>) -> Self {
        let mut buffer = Vec::with_capacity(WIDTH * HEIGHT * 4);
        unsafe {
            buffer.set_len(buffer.capacity());
        }
        buffer.fill(0);
        let pen_color = [255,255,255,255];
        Self { buffer, palette, pen_color }
    }

    pub fn update(&mut self) {
        self.buffer.fill(0);
        self.draw_circle(Vec2::new(100.0, 100.0), 10.0);
        self.draw_line(Vec2::new(20.0, 30.0), Vec2::new(150.0, 20.0));
        return;
        for pixel in self.buffer.as_mut_slice().chunks_exact_mut(4) {
            pixel[0] = (pixel[0] as i32 + fastrand::i32(-2..3)).min(255).max(0) as u8;
            pixel[1] = (pixel[1] as i32 + fastrand::i32(-2..3)).min(255).max(0) as u8;
            pixel[2] = (pixel[2] as i32 + fastrand::i32(-2..3)).min(255).max(0) as u8;
            pixel[3] = (pixel[3] as i32 + fastrand::i32(-2..3)).min(255).max(0) as u8;
        }
    }

    pub fn display(&self) {
        use std::io::Write;
        let file = std::fs::File::options()
            .create(true)
            .read(true)
            .write(true)
            .open("/tmp/imagesink")
            .unwrap();
        let size = 640 * 480 * 4;
        file.set_len(size.try_into().unwrap()).unwrap();
        let mut mmap = unsafe { memmap2::MmapMut::map_mut(&file).unwrap() };
        if let Some(err) = mmap.lock().err() {
            panic!("{err}");
        }
        (&mut mmap[..]).write_all(&self.buffer.as_slice());
    }
    fn random(&mut self) {
        for pixel in self.buffer.as_mut_slice().chunks_exact_mut(4) {
            let change = self.palette[fastrand::usize(0..self.palette.len())];
            pixel[0] = change[0];
            pixel[1] = change[1];
            pixel[2] = change[2];
            pixel[3] = change[3];
        }
    }

    fn draw_line(&mut self, from: Vec2, to: Vec2) {
        let delta = to - from;
        let axis_biggest_distance = (delta.x).abs().max((delta.y).abs()) as usize;
        let normalized = delta.normalize();
        for step in 0..axis_biggest_distance {
            let magnitude = step as f32;
            let x = from.x + normalized.x * magnitude;
            let y = from.y + normalized.y * magnitude;
            self.draw_point(Vec2::new(x, y));
        }
    }

    fn draw_circle(&mut self, pos: Vec2, radius: f32) {
        let left_x = (pos.x - radius) as usize;
        let right_x = (pos.x + radius) as usize;
        let top_y = (pos.y - radius) as usize;
        let bottom_y = (pos.y + radius) as usize;
        for offset_x in left_x..right_x {
            for offset_y in top_y..bottom_y {
                if ((offset_x as f32 - pos.x as f32).powi(2) + (offset_y as f32 - pos.y as f32).powi(2)).sqrt() < radius {
                    self.draw_point(Vec2::new(offset_x as f32, offset_y as f32));
                }
            }
        }
    }

    fn draw_point(&mut self, pos: Vec2) {
        let buffer_idx = self.idx(pos.x as usize, pos.y as usize);
        self.buffer[buffer_idx] = self.pen_color[0];
        self.buffer[buffer_idx + 1] = self.pen_color[1];
        self.buffer[buffer_idx + 2] = self.pen_color[2];
        self.buffer[buffer_idx + 3] = self.pen_color[3];
    }

    fn idx(&self, x: usize, y: usize) -> usize {
        (x + y * WIDTH) * 4
    }
}

fn hex_to_rgb(hex: &&str) -> [u8; 4] {
    let hex = hex.trim_matches('#');
    [
        u8::from_str_radix(&hex[0..2], 16).unwrap(),
        u8::from_str_radix(&hex[2..4], 16).unwrap(),
        u8::from_str_radix(&hex[4..6], 16).unwrap(),
        255,
    ]
}
