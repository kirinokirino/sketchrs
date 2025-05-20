use memmap2::MmapMut;

use std::fs::File;
use std::io::{Read, Write};

const WIDTH: u16 = 640;
const HEIGHT: u16 = 480;

fn main() {
    let mut file = File::options()
        .create(true)
        .read(true)
        .write(true)
        .open("/tmp/imagesink")
        .unwrap();
    let size = usize::from(WIDTH) * usize::from(HEIGHT) * 4;
    file.set_len(size.try_into().unwrap()).unwrap();
    let result = file.write_all(vec![200; size].as_slice());

    let mut mmap = unsafe { MmapMut::map_mut(&file).unwrap() };
    mmap.lock();
    let mut i = 0;
    loop {
        let data = vec![i; size];
        (&mut mmap[..]).write_all(&data);
        i += 10;
        i %= 80;
        std::thread::sleep(std::time::Duration::from_millis(100));
    }
}
