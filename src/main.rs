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
        }
    }

    pub fn update(&mut self) {
        self.acc = Vec2::new(fastrand::f32() * 4. - 2., fastrand::f32() * 4. - 2.);

        self.history.push(self.pos);
        self.pos += self.vel;
        self.vel += self.acc;

        self.stay_on_canvas();

        self.acc = Vec2::new(0.0, 0.0);
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
        canvas.select_color(0);
        for past in &self.history {
            canvas.draw_point(*past);
        }
        let history = self.history.as_slice().windows(3).rev();
        for w in history.clone().take(10) {
            canvas.draw_curve(w[0], w[1], w[2]);
        }
        canvas.select_color(1);
        for w in history.clone().skip(10).take(10) {
            canvas.draw_curve(w[0], w[1], w[2]);
        }
        canvas.select_color(3);
        for w in history.clone().skip(20).take(10) {
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
        let wanderer = Wanderer::new(Vec2::new((WIDTH / 2) as f32, (HEIGHT / 2) as f32));
        let mut wanderers = Vec::new();
        for _ in 0..6 {
            wanderers.push(wanderer.clone());
        }
        wanderers.push(wanderer);
        Self {
            canvas,
            ffmpeg,
            wanderers,
        }
    }

    pub fn run(&mut self) {
        loop {
            self.update();
            self.draw();
            std::thread::sleep(std::time::Duration::from_millis(30));
        }
    }

    fn update(&mut self) {
        self.wanderers.iter_mut().for_each(|w| w.update());
    }

    fn draw(&mut self) {
        //self.canvas.buffer.fill(0);
        self.canvas.dim(1);
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
        //palette.extend([[0, 0, 0, 0]].repeat(1));
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
        buffer.fill(255);
        let pen_color = [255, 255, 255, 255];
        Self {
            buffer,
            palette,
            pen_color,
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
        for pixel in self.buffer.as_mut_slice().chunks_exact_mut(4) {
            let change = self.palette[fastrand::usize(0..self.palette.len())];
            pixel[0] = change[0];
            pixel[1] = change[1];
            pixel[2] = change[2];
            pixel[3] = change[3];
        }
    }

    fn draw_curve(&mut self, start: Vec2, control: Vec2, end: Vec2) {
        let points = start.distance(control) + control.distance(end);
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

    fn draw_point(&mut self, pos: Vec2) {
        let buffer_idx = self.idx(pos.x as usize, pos.y as usize);
        if (buffer_idx + 3) > self.buffer.len() {
            // TODO err?
            return;
        }
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
