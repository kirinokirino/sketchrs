use std::borrow::BorrowMut;

const WIDTH: usize = 640;
const HEIGHT: usize = 480;

fn main() {
    let mut palette: Vec<_> = ["#df8c00", "#d7ac64", "#36241e", "#4a3d35", "#747769"].iter().map(hex_to_rgb).collect();
    palette.extend([[0,0,0,0]].repeat(30));
    let mut canvas = Canvas::new(palette);
    loop {
        canvas.random();
        canvas.display();
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
}

struct Canvas {
    buffer: Vec<u8>,
    palette: Vec<[u8; 4]>,
}

impl Canvas {
    pub fn new(palette: Vec<[u8; 4]>) -> Self {
        let mut buffer = Vec::with_capacity(WIDTH * HEIGHT * 4);
        unsafe {
            buffer.set_len(buffer.capacity());
        }
        buffer.fill(0);
        Self { buffer, palette }
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

    fn draw_circle(&mut self, x: usize, y: usize, radius: f32) {
        todo!();
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
