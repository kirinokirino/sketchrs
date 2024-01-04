use glam::Vec2;

use std::f32::consts::{PI, TAU};
use std::io::Write;
use std::process::{ChildStdin, Command, Stdio};

const RECORD: bool = false;
const WIDTH: usize = 640;
const HEIGHT: usize = 480;
const PALETTE: [&'static str; 5] = ["#ffffff", "#d31055", "#ea6ba3", "#f3b2d4", "#fff0f0"];

fn main() {
    let mut sketch = Sketch::new();
    sketch.run();
}

struct Image {
    origin: Vec2,
    width: usize,
    height: usize,
    pub data: Vec<usize>,
}

impl Image {
    fn new(origin: Vec2, width: usize, height: usize) -> Self {
        Self {
            origin,
            width,
            height,
            data: vec![0; width * height],
        }
    }
}

struct Flower {
    distance: f32,
    angle: f32,
    scale: f32,
    points: Vec<Vec2>,
    image: Option<Image>,
}

impl Flower {
    pub fn new(scale: f32) -> Self {
        Self {
            distance: 0.0,
            angle: 0.0,
            scale,
            points: Vec::new(),
            image: None,
        }
    }

    pub fn add_point(&mut self) {
        let center = Vec2::new((WIDTH / 2) as f32, (HEIGHT / 2) as f32);

        let pos = center + Vec2::from_angle(self.angle.to_radians()) * self.distance;
        self.points.push(pos);

        self.angle += 137.5;
        if self.angle >= 360.0 {
            self.angle -= 360.0;
        }
        self.distance += 0.2;
    }

    pub fn draw(&self, canvas: &mut Canvas) {
        if let Some(image) = &self.image {
            println!("IMAGE {} {} {}", image.origin, image.width, image.height);
            for (idx, palette_color) in image.data.iter().enumerate() {
                let x = image.origin.x + (idx % image.width) as f32;
                let y = image.origin.y + (idx / image.width) as f32;
                let pos = Vec2::new(x, y);
                if *palette_color == 0 {
                    canvas.transparent();
                } else {
                    canvas.select_color(*palette_color as u8);
                }
                canvas.draw_point(pos)
            }
        }
    }

    fn generate_image(&mut self) {
        let mut min_x = 9999.0;
        let mut min_y = 9999.0;
        let mut max_x = -9999.0;
        let mut max_y = -9999.0;
        for point in &self.points {
            let (x, y) = (point.x, point.y);
            if x < min_x {
                min_x = x;
            }
            if x > max_x {
                max_x = x;
            }
            if y < min_y {
                min_y = y;
            }
            if y > max_y {
                max_y = y;
            }
        }

        let origin = Vec2::new(min_x, min_y);
        let width = (max_x - min_x).ceil().max(1.0) as usize;
        let height = (max_y - min_y).ceil().max(1.0) as usize;
        let mut image = Image::new(origin, width, height);

        let voronoi = true;
        let center = Vec2::new((width / 2) as f32, (height / 2) as f32);
        if voronoi {
            for x in 0..width {
                for y in 0..height {
                    let point = Vec2::new(x as f32, y as f32);
                    if point.distance(center) > self.distance {
                        continue;
                    }
                    let closest = find_closest(point + origin, &self.points);
                    image.data[x + y * width] = 1 + closest % 4;
                }
            }
        } else {
            for (i, point) in self.points.iter().enumerate() {
                let color = i % 4;
                let local_point = *point - origin;
                let idx = local_point.x as usize + local_point.y as usize * width;
                image.data[idx] = 1 + color;
            }
        }
        self.image = Some(image);
    }
}

fn find_closest(point: Vec2, points: &[Vec2]) -> usize {
    let mut min_distance = 99999.0;
    let mut idx = 0;
    for (i, other) in points.iter().enumerate() {
        let dist = point.distance(*other);
        if dist < min_distance {
            min_distance = dist;
            idx = i
        }
    }
    idx
}

struct Sketch {
    canvas: Canvas,
    ffmpeg: Option<ChildStdin>,
    cycle: usize,
    flower: Flower,
}

impl Sketch {
    pub fn new() -> Self {
        let ffmpeg = Self::ffmpeg();
        let canvas = Self::canvas();
        let flower = Flower::new(2.0);

        Self {
            canvas,
            ffmpeg,
            cycle: 0,
            flower,
        }
    }

    pub fn run(&mut self) {
        loop {
            self.update();
            self.draw();
            //std::thread::sleep(std::time::Duration::from_millis(30));
            self.cycle += 1;
        }
    }

    fn update(&mut self) {
        if self.cycle % 1 == 0 {
            self.flower.add_point();
            self.flower.add_point();
            self.flower.add_point();
            self.flower.add_point();
            self.flower.add_point();
            self.flower.add_point();
            self.flower.add_point();
            self.flower.add_point();
            self.flower.add_point();
            self.flower.add_point();
            self.flower.generate_image();
        }
    }

    fn draw(&mut self) {
        // self.canvas.blend_mode = BlendMode::Blend;
        // self.canvas.pen_color = hex_to_rgb(&PALETTE[PALETTE.len() - 1]);
        // self.canvas.pen_color[3] = 20;
        // self.canvas.draw_square(
        //     Vec2::new(0.0, 0.0),
        //     Vec2::new((WIDTH) as f32, HEIGHT as f32),
        // );
        // self.canvas.blend_mode = BlendMode::Replace;
        self.canvas.buffer.fill(0);
        // self.canvas.dim(10);
        //self.canvas.random();
        self.canvas.pen_color = hex_to_rgb(&PALETTE[0]);
        self.flower.draw(&mut self.canvas);

        if RECORD {
            self.ffmpeg
                .as_mut()
                .map(|ffmpeg| ffmpeg.write_all(&self.canvas.buffer.as_slice()));
        }
        self.canvas.display();
    }

    fn canvas() -> Canvas {
        let mut palette: Vec<_> = PALETTE.iter().map(hex_to_rgb).collect();
        palette.extend([[0, 0, 0, 0]].repeat(1));
        Canvas::new(palette)
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

enum BlendMode {
    Replace,
    Blend,
}

struct Canvas {
    pub buffer: Vec<u8>,
    palette: Vec<[u8; 4]>,
    pub pen_color: [u8; 4],
    blend_mode: BlendMode,
}

impl Canvas {
    pub fn new(palette: Vec<[u8; 4]>) -> Self {
        let mut buffer = vec![255u8; WIDTH * HEIGHT * 4];
        let pen_color = [255, 255, 255, 255];
        Self {
            buffer,
            palette,
            pen_color,
            blend_mode: BlendMode::Replace,
        }
    }

    pub fn select_color(&mut self, color: u8) {
        self.pen_color = self.palette[color as usize % self.palette.len()]
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

    fn random(&mut self) {
        for i in 0..self.buffer.len() / 4 {
            let mut change = self.palette[fastrand::usize(0..self.palette.len())];
            change[3] = (change[3] as f32 * 0.05) as u8;
            self.pen_color = change;
            self.point_blend(i * 4);
        }
    }

    pub fn draw_buffer(&mut self, image: &Image) {
        todo!();
    }

    pub fn fluffy_line(&mut self, start: Vec2, end: Vec2) {
        let middle = (start + end) / 2.0;
        let randomness = 12.0;
        for _ in 0..4 {
            self.select_color(fastrand::u8(0..PALETTE.len() as u8));
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
                if ((offset_x as f32 - pos.x as f32).powi(2)
                    + (offset_y as f32 - pos.y as f32).powi(2))
                .sqrt()
                    < radius
                {
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

    fn draw_square(&mut self, top_left: Vec2, bottom_right: Vec2) {
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
