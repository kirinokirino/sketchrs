## A window that's displaying the contents of `/tmp/imagesink` file

Run this app in tandem with ![sketchrs](https://github.com/kirinokirino/sketchrs) or some other program that dumps its visual buffer in `/tmp/imagesink`, and you'll get a graphical application witout handling window creation.

Only works on linux, there is probably a windows feature that I can enable.

To install, run `cargo install --path=.`. Then you can run it with your app, for example `cargo run & imagesink`
