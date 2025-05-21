use glam::Vec2;

use std::f32::consts::PI;
use std::io::Write;
use std::process::{ChildStdin, Command, Stdio};

const RECORD: bool = false;
const WIDTH: usize = 640;
const HEIGHT: usize = 480;

fn main() {
    let mut sketch = Sketch::new();
    sketch.run();
}

struct Sketch {
    canvas: Canvas,
    ffmpeg: Option<ChildStdin>,
    palette: Vec<[u8; 4]>,
    dispersion_matrix: Vec<f32>,
}

impl Sketch {
    pub fn new() -> Self {
        let palette = [
            "#0e0e12", "#1a1a24", "#333346", "#535373", "#8080a4", "#a6a6bf", "#c1c1d2", "#e6e6ec",
        ]
        .iter()
        .map(hex_to_rgb)
        .collect::<Vec<_>>();

        let dispersion_matrix = vec![0.11, 0.77, 0.44, 0.55, 0.88, 0.33, 0.66, 0.22, 1.0];

        let ffmpeg = Self::ffmpeg();
        let canvas = Self::canvas();

        Self {
            canvas,
            ffmpeg,
            palette,
            dispersion_matrix,
        }
    }

    pub fn run(&mut self) {
        let target = Vec2::new(WIDTH as f32 / 2.0, HEIGHT as f32 * 2.0);
        let min_distance = Vec2::new(WIDTH as f32 / 2.0, HEIGHT as f32).distance_squared(target);
        let max_distance = Vec2::new(0.0, 0.0).distance_squared(target);
        for x in 0..WIDTH {
            for y in 0..HEIGHT {
                let pos = Vec2::new(x as f32, y as f32);
                let distance = pos.distance_squared(target);
                let gradient = map(distance, min_distance, max_distance, 0.0, 1.0);
                let color = self.dispersed_pixel_color(pos, gradient);
                self.canvas.pen_color = color;
                self.canvas.draw_point(pos);
            }
        }

        loop {
            self.draw();
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
    }

    fn draw(&mut self) {
        if RECORD {
            self.ffmpeg
                .as_mut()
                .map(|ffmpeg| ffmpeg.write_all(&self.canvas.buffer.as_slice()));
        }
        self.canvas.display();
    }

    fn dispersed_pixel_color(&self, pos: Vec2, brightness: f32) -> [u8; 4] {
        let x = pos.x as usize;
        let y = pos.y as usize;
        let index = (x.overflowing_sub(y * 3).0) % self.dispersion_matrix.len();
        let dispersion = self.dispersion_matrix[index] * self.brightness_step();
        let brightness = (brightness - self.brightness_step()) + dispersion;
        self.color_by_brightness(brightness)
    }

    fn color_by_brightness(&self, brightness: f32) -> [u8; 4] {
        let step = self.brightness_step();
        let index = (brightness / step).round() as usize;
        if index >= self.palette.len() {
            return self.palette.last().unwrap().clone();
        }
        self.palette[index]
    }

    fn brightness_step(&self) -> f32 {
        let step = 1.0 / self.palette.len() as f32;
        step
    }

    fn canvas() -> Canvas {
        Canvas::new()
    }

    fn ffmpeg() -> Option<ChildStdin> {
        let ffmpeg_command = "/usr/bin/ffmpeg";
        let args = "-y -f rawvideo -vcodec rawvideo -s 640x480 -pix_fmt rgba -r 30 -i - -an -vcodec h264 -pix_fmt yuv420p -crf 15 /home/k/Media/video.mp4".split(' ');
        if RECORD {
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
        }
    }
}
#[derive(PartialEq)]
enum BlendMode {
    Replace,
    Blend,
}

struct Canvas {
    pub buffer: Vec<u8>,
    pub pen_color: [u8; 4],
    blend_mode: BlendMode,
}

impl Canvas {
    pub fn new() -> Self {
        let buffer = vec![255u8; WIDTH * HEIGHT * 4];
        let pen_color = [255, 255, 255, 255];
        Self {
            buffer,
            pen_color,
            blend_mode: BlendMode::Replace,
        }
    }

    pub fn transparent(&mut self) {
        self.pen_color = [0, 0, 0, 0];
    }

    pub fn dim(&mut self, value: i16) {
        self.buffer.iter_mut().for_each(|v| {
            let new = (*v as i16 + value).min(255).max(0) as u8;
            *v = new;
        });
    }

    pub fn display(&self) {
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
        let _ = (&mut mmap[..]).write_all(&self.buffer.as_slice());
    }

    pub fn fluffy_line(&mut self, start: Vec2, end: Vec2) {
        let middle = (start + end) / 2.0;
        let randomness = 12.0;
        for _ in 0..4 {
            self.draw_curve(
                start
                    + Vec2::new(
                        (fastrand::f32() - 0.5) * randomness,
                        (fastrand::f32() - 0.5) * randomness,
                    ),
                middle
                    + Vec2::new(
                        (fastrand::f32() - 0.5) * randomness,
                        (fastrand::f32() - 0.5) * randomness,
                    ),
                end + Vec2::new(
                    (fastrand::f32() - 0.5) * randomness / 2.0,
                    (fastrand::f32() - 0.5) * randomness / 2.0,
                ),
            )
        }
    }

    fn draw_curve(&mut self, start: Vec2, control: Vec2, end: Vec2) {
        let points = start.distance(control) + control.distance(end) + end.distance(start);
        for i in 1..points as usize {
            let proportion = i as f32 / points;
            let path1 = control - start;
            let point1 = start + path1 * proportion;
            let path2 = end - control;
            let point2 = control + path2 * proportion;
            let path3 = point2 - point1;
            let point3 = point1 + path3 * proportion;
            self.draw_point(point3);
        }
    }

    fn draw_curve_dots(&mut self, start: Vec2, control: Vec2, end: Vec2, dots_density: f32) {
        let points =
            (start.distance(control) + control.distance(end) + end.distance(start)) * dots_density;
        for i in 1..points as usize {
            let proportion = i as f32 / points;
            let path1 = control - start;
            let point1 = start + path1 * proportion;
            let path2 = end - control;
            let point2 = control + path2 * proportion;
            let path3 = point2 - point1;
            let point3 = point1 + path3 * proportion;
            self.draw_point(point3);
        }
    }

    fn draw_line(&mut self, from: Vec2, to: Vec2) {
        let delta = to - from;
        let normalized = delta.normalize();
        for step in 0..delta.length() as usize {
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
        for offset_x in left_x..=right_x {
            for offset_y in top_y..=bottom_y {
                let distance = (offset_x as f32 - pos.x as f32).powi(2)
                    + (offset_y as f32 - pos.y as f32).powi(2);
                if distance.sqrt() < radius {
                    let brightness = map(distance.sqrt(), 0.0, radius, 0.0, 1.0);
                    self.draw_point(Vec2::new(offset_x as f32, offset_y as f32));
                }
            }
        }
    }

    fn draw_arc(&mut self, center: Vec2, radius: f32, angle_spread: f32, direction: f32) {
        let steps = (2.0 * radius * PI) as usize + 1;
        let start = direction - angle_spread / 2.0;
        let end = direction + angle_spread / 2.0;
        for step in 0..steps {
            let arc_point =
                Vec2::from_angle(map(step as f32, 0.0, steps as f32, start, end)) * radius + center;
            self.draw_point(arc_point);
        }
    }

    fn draw_rectangle(&mut self, top_left: Vec2, bottom_right: Vec2) {
        for offset_x in top_left.x as usize..=bottom_right.x as usize {
            for offset_y in top_left.y as usize..=bottom_right.y as usize {
                self.draw_point(Vec2::new(offset_x as f32, offset_y as f32));
            }
        }
    }

    fn draw_point(&mut self, pos: Vec2) {
        if pos.x >= 640.0 || pos.x < 0.0 || pos.y >= 480.0 || pos.y < 0.0 {
            return;
        }
        let buffer_idx = self.idx(pos.x as usize, pos.y as usize);
        // if (buffer_idx + 3) > self.buffer.len() {
        //     // TODO err?
        //     return;
        // }
        match self.blend_mode {
            BlendMode::Replace => self.point_replace(buffer_idx),
            BlendMode::Blend => self.point_blend(buffer_idx),
        }
    }

    fn point_blend(&mut self, buffer_idx: usize) {
        let [r, g, b, a] = self.pen_color;

        if a == 0 {
            return;
        } else if a == 255 {
            self.point_replace(buffer_idx);
            return;
        }

        let mix = a as f32 / 255.0;
        let [dst_r, dst_g, dst_b, dst_a] = [
            self.buffer[buffer_idx] as f32,
            self.buffer[buffer_idx + 1] as f32,
            self.buffer[buffer_idx + 2] as f32,
            self.buffer[buffer_idx + 3] as f32,
        ];

        self.buffer[buffer_idx] = ((r as f32 * mix) + (dst_r * (1.0 - mix))) as u8;
        self.buffer[buffer_idx + 1] = ((g as f32 * mix) + (dst_g * (1.0 - mix))) as u8;
        self.buffer[buffer_idx + 2] = ((b as f32 * mix) + (dst_b * (1.0 - mix))) as u8;
        self.buffer[buffer_idx + 3] = ((a as f32 * mix) + (dst_a * (1.0 - mix))) as u8;
    }

    fn point_replace(&mut self, buffer_idx: usize) {
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

pub fn map(value: f32, start1: f32, stop1: f32, start2: f32, stop2: f32) -> f32 {
    (value - start1) / (stop1 - start1) * (stop2 - start2) + start2
}

pub fn lerp(t: f32, start: f32, stop: f32) -> f32 {
    start + t * (stop - start)
}
