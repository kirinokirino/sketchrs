use std::f32::consts::TAU;
use std::io::Write;
use std::process::{ChildStdin, Command, Stdio};

const RECORD: bool = true;
const WIDTH: usize = 640;
const HEIGHT: usize = 480;
const PALETTE: [&'static str; 10] = [
    "#c76e86", "#9c2940", "#761831", "#420424", "#090b19", "#091633", "#14213e", "#24454c",
    "#316264", "#407b70",
];

use glam::Vec2;

fn main() {
    let mut sketch = Sketch::new();
    sketch.run();
}
struct Cloud {
    points: Vec<Vec2>,
}

impl Cloud {
    pub fn new(pos: Vec2, size: Vec2) -> Self {
        let randomness = 0.5;
        let max_points = 1000;
        let angle_step = TAU / max_points as f32;
        let mut points = Vec::with_capacity(max_points);
        for i in 0..max_points {
            let rotation = Vec2::from_angle(angle_step * i as f32);
            let point = Vec2::ONE.rotate(rotation);
            let point_scaled = point * size + pos;
            let point_randomized = point_scaled
                + point_scaled
                    * Vec2::new(fastrand::f32() * randomness, fastrand::f32() * randomness);
            points.push(point_randomized);
        }
        Self { points }
    }

    pub fn draw(&self, offset: u32, canvas: &mut Canvas) {
        let offset_max = 30;
        for (i, w) in self.points.as_slice().windows(6 + offset_max).enumerate() {
            let offset = (offset as usize % offset_max) as usize;
            if (i % offset_max) != offset {
                continue;
            }
            let color = ((self.points[0].y / 480.0) * PALETTE.len() as f32) as u8;
            canvas.select_color(color);
            canvas.draw_curve(w[0 + offset], w[2 + offset], w[4 + offset]);
        }
    }

    pub fn update(&mut self) {
        for point in &mut self.points {
            *point += Vec2::new((fastrand::f32() - 0.5) * 6.0, (fastrand::f32() - 0.5) * 6.0);
        }
    }
}

struct Sketch {
    canvas: Canvas,
    ffmpeg: Option<ChildStdin>,

    clouds: Vec<Cloud>,
    cycle: u32,
}

impl Sketch {
    pub fn new() -> Self {
        let ffmpeg = Self::ffmpeg();
        let canvas = Self::canvas();
        let mut clouds = Vec::new();
        let size: f32 = 40.0;
        for _ in 0..40 {
            let pos = Vec2::new(fastrand::f32() * 640.0, -fastrand::f32() * 240.0 + 240.0);
            let size = Vec2::new(
                fastrand::f32() * size + size * 0.3,
                fastrand::f32() / 3.0 * size + size * 0.3,
            );
            clouds.push(Cloud::new(pos, size));
        }
        Self {
            canvas,
            ffmpeg,
            clouds,
            cycle: 0,
        }
    }

    pub fn run(&mut self) {
        loop {
            self.update();
            self.draw();
            std::thread::sleep(std::time::Duration::from_millis(30));
            self.cycle += 1;
        }
    }

    fn update(&mut self) {
        self.clouds.iter_mut().for_each(|cloud| cloud.update());
    }

    fn draw(&mut self) {
        self.canvas.blend_mode = BlendMode::Blend;
        self.canvas.pen_color = hex_to_rgb(&PALETTE[PALETTE.len() - 1]);
        self.canvas.pen_color[3] = 20;
        self.canvas.draw_square(
            Vec2::new(0.0, 0.0),
            Vec2::new((WIDTH) as f32, HEIGHT as f32),
        );
        self.canvas.blend_mode = BlendMode::Replace;
        // self.canvas.buffer.fill(0);
        self.canvas.dim(10);
        //self.canvas.random();
        self.clouds
            .iter()
            .for_each(|cloud| cloud.draw(self.cycle, &mut self.canvas));

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
    pub pen_color: [u8; 4],
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
