//! File-Based AEC Processing Tool
//!
//! This example acts as a command-line utility to perform acoustic echo cancellation on WAV files.
//! It takes a far-end signal file and a microphone signal file as input, processes them,
//! and saves the echo-cancelled signal to a new WAV file.
//!
//! This is the most "professional" and practical example, as it allows you to test the AEC
//! with your own audio recordings.
//!
//! ## How to Run
//!
//! 1. Place your far-end and microphone WAV files in the root of the `fdaf_aec` directory.
//!    (The files should be single-channel, 16kHz, 16-bit PCM for best results).
//! 2. Run the command below, replacing the file names with your own.
//!
//! ```sh
//! cargo run --example file_based_aec --release -- \
//!   --farend your_farend_file.wav \
//!   --mic your_mic_file.wav \
//!   --output your_output_file.wav
//! ```
//! The `--release` flag is recommended for faster processing.

use fdaf_aec::FdafAec;
use hound::{WavReader, WavWriter, WavSpec};
use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Path to the far-end signal WAV file (e.g., loudspeaker audio).
    #[clap(long, value_parser)]
    farend: PathBuf,

    /// Path to the microphone signal WAV file (contains echo and near-end speech).
    #[clap(long, value_parser)]
    mic: PathBuf,

    /// Path to save the processed output WAV file.
    #[clap(long, value_parser)]
    output: PathBuf,

    /// Step size (learning rate) for the adaptive filter.
    #[clap(long, value_parser, default_value_t = 0.02)]
    step_size: f32,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    println!("--- Running File-Based AEC ---");
    println!("- Far-end file: {}", args.farend.display());
    println!("- Mic file:     {}", args.mic.display());
    println!("- Output file:  {}", args.output.display());
    println!("- Step size:    {}", args.step_size);

    const FFT_SIZE: usize = 1024;
    const FRAME_SIZE: usize = FFT_SIZE / 2;

    // --- 1. Load WAV files ---
    let (far_signal, spec) = read_wav(&args.farend)?;
    let (mic_signal, _) = read_wav(&args.mic)?;
    
    if spec.channels != 1 || spec.sample_rate != 16000 {
        eprintln!("Warning: For best results, input WAV files should be mono, 16kHz.");
        eprintln!("Current spec: {} channels, {} Hz", spec.channels, spec.sample_rate);
    }
    
    // --- 2. Initialize AEC ---
    let mut aec = FdafAec::new(FFT_SIZE, args.step_size);

    // --- 3. Process Signals Frame by Frame ---
    let mut processed_signal = Vec::new();
    let num_samples = mic_signal.len().min(far_signal.len());

    for i in (0..num_samples).step_by(FRAME_SIZE) {
        if i + FRAME_SIZE > num_samples { break; }

        let far_frame = &far_signal[i..i + FRAME_SIZE];
        let mic_frame = &mic_signal[i..i + FRAME_SIZE];

        let output_frame = aec.process(far_frame, mic_frame);
        processed_signal.extend_from_slice(&output_frame);
    }
    
    // --- 4. Save Output WAV File ---
    let mut writer = WavWriter::create(&args.output, spec)?;
    for &sample in processed_signal.iter() {
        writer.write_sample((sample.clamp(-1.0, 1.0) * i16::MAX as f32) as i16)?;
    }
    writer.finalize()?;

    println!("\nProcessing complete!");
    println!("Output saved to '{}'", args.output.display());
    Ok(())
}

/// Reads a WAV file and converts its samples to a Vec<f32> in the range [-1.0, 1.0].
fn read_wav(path: &PathBuf) -> Result<(Vec<f32>, WavSpec), Box<dyn std::error::Error>> {
    let mut reader = WavReader::open(path)?;
    let spec = reader.spec();
    let max_val = 2_i32.pow(spec.bits_per_sample as u32 - 1) as f32;

    let samples = reader.samples::<i32>()
        .map(|s| s.unwrap() as f32 / max_val)
        .collect();
    
    Ok((samples, spec))
} 