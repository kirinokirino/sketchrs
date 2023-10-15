use image::io::Reader as ImageReader;
use image::{Pixel, Rgba, RgbaImage};
use std::f32::consts::TAU;
use std::io::Write;
use std::process::{ChildStdin, Command, Stdio};

const RECORD: bool = false;
const WIDTH: usize = 640;
const HEIGHT: usize = 480;
const PALETTE: [&'static str; 32] = [
    "#6074ab", "#6b9acf", "#8bbde6", "#aae0f3", "#c8eded", "#faffe0", "#dde6e0", "#b4bec2",
    "#949da8", "#7a7a99", "#5b5280", "#4e3161", "#421e42", "#612447", "#7a3757", "#96485b",
    "#bd6868", "#d18b79", "#dbac8c", "#e6cfa1", "#e7ebbc", "#b2dba0", "#87c293", "#70a18f",
    "#637c8f", "#b56e75", "#c98f8f", "#dfb6ae", "#edd5ca", "#bd7182", "#9e5476", "#753c6a",
];

const SCALE: usize = 40;
const CELLS_WIDTH: usize = WIDTH / SCALE;
const CELLS_HEIGHT: usize = HEIGHT / SCALE;

use glam::Vec2;

fn main() {
    let mut sketch = Sketch::new();
    sketch.run();
}

#[derive(Clone)]
struct Cell {
    pos: Vec2,
    color: [u8; 4],
}

impl Cell {
    pub fn new(pos: Vec2, color: [u8; 4]) -> Self {
        Self { pos, color }
    }

    pub fn draw(&self, offset: u32, canvas: &mut Canvas) {
        //
    }

    pub fn update(&mut self) {
        self.pos += Vec2::new(fastrand::f32() - 0.5, fastrand::f32() - 0.5);
    }
}

struct Cells {
    grid: [[Vec<Cell>; CELLS_WIDTH]; CELLS_HEIGHT],
    len: usize,
}

impl Cells {
    pub fn new() -> Self {
        let mut grid: [[std::mem::MaybeUninit<Vec<Cell>>; CELLS_WIDTH]; CELLS_HEIGHT] =
            unsafe { std::mem::MaybeUninit::uninit().assume_init() };
        for y in 0..CELLS_HEIGHT {
            for x in 0..CELLS_WIDTH {
                grid[y][x] = std::mem::MaybeUninit::new(vec![Cell::new(
                    Vec2::new(
                        x as f32 * (SCALE as f32) + (SCALE as f32) / 2.0,
                        y as f32 * (SCALE as f32) + (SCALE as f32) / 2.0,
                    ),
                    [0, 0, 0, 0], //[0, 181, 226, 255],
                )]);
            }
        }
        let len = CELLS_WIDTH * CELLS_HEIGHT;
        Self {
            grid: unsafe {
                std::mem::transmute::<_, [[Vec<Cell>; CELLS_WIDTH]; CELLS_HEIGHT]>(grid)
            },
            len,
        }
    }

    pub fn push(&mut self, cell: Cell) {
        let pos = cell.pos;
        let x = (pos.x / (SCALE as f32)) as usize;
        let y = (pos.y / (SCALE as f32)) as usize;
        let other = self.closest(pos);
        if other.pos.distance(pos) < 1.5 {
            self.closest_brighten(pos, cell.color);
        } else {
            self.grid[y][x].push(cell);
            self.len += 1;
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn closest_brighten(&mut self, pos: Vec2, color: [u8; 4]) {
        let cells =
            &mut self.grid[(pos.y / (SCALE as f32)) as usize][(pos.x / (SCALE as f32)) as usize];
        let mut closest = &mut cells[0].clone();
        let mut closest_distance = 9999f32;

        for (_i, cell) in cells.iter_mut().enumerate() {
            let dist = cell.pos.distance(pos);
            if dist < closest_distance {
                closest_distance = dist;
                closest = cell;
            }
        }
        let c = &mut closest.color;
        c[0] = (c[0].max(color[0]) as u16 + 1) as u8;
        c[1] = (c[1].max(color[1]) as u16 + 1) as u8;
        c[2] = (c[2].max(color[2]) as u16 + 1) as u8;
    }

    pub fn closest(&self, pos: Vec2) -> &Cell {
        let cells =
            &self.grid[(pos.y / (SCALE as f32)) as usize][(pos.x / (SCALE as f32)) as usize];
        let mut closest = &cells[0];
        let mut closest_distance = 9999f32;

        for cell in cells.iter() {
            let dist = cell.pos.distance(pos);
            if dist < closest_distance {
                closest_distance = dist;
                closest = cell;
            }
        }
        closest
    }
}

struct Sketch {
    canvas: Canvas,
    ffmpeg: Option<ChildStdin>,

    cells: Cells,
    castle: Vec<(Vec2, u8)>,
    cycle: u32,
    max_brightness: u8,
}

impl Sketch {
    pub fn new() -> Self {
        let ffmpeg = Self::ffmpeg();
        let canvas = Self::canvas();
        let mut cells = Cells::new();
        let img: image::ImageBuffer<Rgba<u8>, Vec<u8>> = ImageReader::open("assets/castle.png")
            .unwrap()
            .decode()
            .unwrap()
            .into_rgba8();
        let mut castle = Vec::new();
        let mut max_brightness = 0;
        for (x, y, pixel) in img.enumerate_pixels() {
            if pixel.0[3] > 10 {
                let c = pixel.channels();
                let unred = 255 - c[0];
                castle.push((Vec2::new(x as f32, y as f32), unred));
                if unred > max_brightness {
                    max_brightness = unred;
                }
            }
        }

        Self {
            canvas,
            ffmpeg,
            cells,
            castle,
            cycle: 0,
            max_brightness,
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
        //self.cells.iter_mut().for_each(|cloud| cloud.update());
        if self.cells.len() > 100000 {
            return;
        }
        for i in 0..32 * 4 {
            let (pos, brightness) = self.castle[fastrand::usize(..self.castle.len())];
            let color = self.canvas.palette[i / 4];
            let dim: f32 = brightness as f32 / self.max_brightness as f32;
            let new_color = [
                (color[0] as f32 * dim) as u8,
                (color[1] as f32 * dim) as u8,
                (color[2] as f32 * dim) as u8,
                color[3],
            ];
            self.cells.push(Cell::new(pos, new_color));
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
        // self.canvas.buffer.fill(0);
        // self.canvas.dim(10);
        //self.canvas.random();
        for (i, pixel) in self
            .canvas
            .buffer
            .as_mut_slice()
            .chunks_exact_mut(4)
            .enumerate()
        {
            let x = i % 640;
            let y = i / 640;
            let pix_pos = Vec2::new(x as f32, y as f32);

            let closest = self.cells.closest(pix_pos);
            let dist = closest.pos.distance(pix_pos);
            let c = closest.color;
            pixel[0] = c[0]; //(self.cells[closest].color[0] as f32 * (closest_distance / 50.0)) as u8;
            pixel[1] = c[1]; //(self.cells[closest].color[1] as f32 * (closest_distance / 50.0)) as u8;
            pixel[2] = c[2]; //(self.cells[closest].color[2] as f32 * (closest_distance / 50.0)) as u8;
            pixel[3] = c[3]; //((SCALE as f32 - dist) * (255.0 / SCALE as f32)) as u8;
        }

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
