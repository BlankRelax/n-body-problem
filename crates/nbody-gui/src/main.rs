//! Native entry point. The application lives in the library so that the same
//! code runs in a desktop window and, compiled to WebAssembly, in a web page
//! (where the library's exported `start` is the entry point instead).
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result {
    nbody_gui::run_native()
}

#[cfg(target_arch = "wasm32")]
fn main() {}
