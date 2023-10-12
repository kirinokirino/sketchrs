use std::io::Write;
use std::process::{ChildStdin, Command, Stdio};

const RECORD: bool = false;
const WIDTH: usize = 640;
const HEIGHT: usize = 480;

use glam::Vec2;

fn main() {
    let mut sketch = Sketch::new();
    sketch.run();
}

#[derive(Debug, Clone)]
struct Wanderer {
    pos: Vec2,
    vel: Vec2,
    acc: Vec2,

    history: Vec<Vec2>,
    lifetime: u32,
}

impl Wanderer {
    pub fn new(pos: Vec2) -> Self {
        let vel = Vec2::new(0.0, 0.0);
        let acc = Vec2::new(0.0, 0.0);
        let history = Vec::new();
        Self {
            pos,
            vel,
            acc,
            history,
            lifetime: 0,
        }
    }

    pub fn reset(&mut self) {
        self.pos = Vec2::new(320.0, 240.0);
        self.vel = Vec2::ZERO;
        self.acc = Vec2::ZERO;
        self.history.clear();
    }

    pub fn update(&mut self) {
        self.acc = Vec2::new(fastrand::f32() * 1. - 0.5, fastrand::f32() * 1. - 0.52);

        if self.lifetime % 1 == 0 {
            self.history.push(self.pos);
        }
        self.pos += self.vel;
        self.vel += self.acc;

        //self.stay_on_canvas();
        if self.pos.distance(Vec2::new(320.0, 240.0)) > 640.0 {
            self.reset();
        }

        self.acc = Vec2::new(0.0, 0.0);
        self.lifetime += 1;
    }

    fn stay_on_canvas(&mut self) {
        if self.pos.x <= 0.0 {
            self.pos.x = WIDTH as f32;
        } else if self.pos.x >= WIDTH as f32 {
            self.pos.x = 0.0;
        }
        if self.pos.y <= 0.0 {
            self.pos.y = HEIGHT as f32;
        } else if self.pos.y >= HEIGHT as f32 {
            self.pos.y = 0.0;
        }
    }

    pub fn draw(&self, canvas: &mut Canvas) {
        // for past in &self.history {
        //     canvas.draw_point(*past);
        // }
        let history = self.history.as_slice().windows(3).rev();
        for (i, w) in history.clone().take(10).enumerate() {
            canvas.select_color((i % 5) as u8);
            canvas.draw_curve(w[0], w[1], w[2]);
        }
    }
}

struct Sketch {
    canvas: Canvas,
    ffmpeg: Option<ChildStdin>,

    wanderers: Vec<Wanderer>,
}

impl Sketch {
    pub fn new() -> Self {
        let ffmpeg = Self::ffmpeg();
        let canvas = Self::canvas();
        Self {
            canvas,
            ffmpeg,
            wanderers: Vec::new(),
        }
    }

    pub fn run(&mut self) {
        let mut cycle = 0;
        let wanderer = Wanderer::new(Vec2::new((WIDTH / 2) as f32, (HEIGHT / 2) as f32));

        loop {
            if cycle < 100 {
                self.wanderers.push(wanderer.clone());
            }
            self.update();
            self.draw();
            std::thread::sleep(std::time::Duration::from_millis(30));
            cycle += 1;
        }
    }

    fn update(&mut self) {
        self.wanderers.iter_mut().for_each(|w| w.update());
    }

    fn draw(&mut self) {
        //self.canvas.select_color(1);
        // self.canvas.draw_square(
        //     Vec2::new(0.0, 0.0),
        //     Vec2::new((WIDTH / 2) as f32, HEIGHT as f32),
        // );
        self.canvas.buffer.fill(0);
        //self.canvas.dim(10);
        //self.canvas.random();
        self.wanderers.iter().for_each(|w| w.draw(&mut self.canvas));

        if RECORD {
            self.ffmpeg
                .as_mut()
                .map(|ffmpeg| ffmpeg.write_all(&self.canvas.buffer.as_slice()));
        }
        self.canvas.display();
    }

    fn canvas() -> Canvas {
        let mut palette: Vec<_> = ["#df8c00", "#d7ac64", "#36241e", "#4a3d35", "#747769"]
            .iter()
            .map(hex_to_rgb)
            .collect();
        palette.extend([[0, 0, 0, 0]].repeat(1));
        Canvas::new(palette)
    }

    fn ffmpeg() -> Option<ChildStdin> {
        let ffmpeg_command = "/usr/bin/ffmpeg";
        let args = "-y -f rawvideo -vcodec rawvideo -s 640x480 -pix_fmt rgba -r 30 -i - -an -vcodec h264 -pix_fmt yuv420p -crf 15 /home/kirinokirino/Media/video.mp4".split(' ');
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
    pen_color: [u8; 4],
    blend_mode: BlendMode,
}

impl Canvas {
    pub fn new(palette: Vec<[u8; 4]>) -> Self {
        let mut buffer = Vec::with_capacity(WIDTH * HEIGHT * 4);
        unsafe {
            buffer.set_len(buffer.capacity());
        }
        buffer.fill(255);
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
        self.buffer[buffer_idx] = ((g as f32 * mix) + (dst_g * (1.0 - mix))) as u8;
        self.buffer[buffer_idx] = ((b as f32 * mix) + (dst_b * (1.0 - mix))) as u8;
        self.buffer[buffer_idx] = ((a as f32 * mix) + (dst_a * (1.0 - mix))) as u8;
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
