//! AEC Simulation with Generated Signals
//!
//! This example creates a more realistic AEC test scenario compared to `basic_simulation`.
//! It generates white noise for the far-end signal and a sine wave for the near-end signal.
//! A simple Finite Impulse Response (FIR) filter is used to simulate a basic
//! Room Impulse Response (RIR), creating a more realistic echo.
//! The script saves the far-end, near-end, microphone (before AEC), and output (after AEC)
//! signals as WAV files for auditory comparison.
//!
//! ## How to Run
//!
//! ```sh
//! cargo run --example generated_signal_aec --release
//! ```
//! The `--release` flag is recommended for faster processing.

use fdaf_aec::FdafAec;
use hound;
use rand::{Rng, thread_rng};

const SAMPLE_RATE: u32 = 16000;
const DURATION_S: u32 = 5;

fn main() {
    println!("--- Running AEC Simulation with Generated Signals ---");

    const FFT_SIZE: usize = 1024;
    const FRAME_SIZE: usize = FFT_SIZE / 2;
    const STEP_SIZE: f32 = 0.02;
    const TOTAL_SAMPLES: usize = (SAMPLE_RATE * DURATION_S) as usize;

    let mut aec = FdafAec::new(FFT_SIZE, STEP_SIZE);

    // --- 1. Signal Generation ---
    let mut rng = thread_rng();

    // Far-end signal (white noise) to simulate a generic, broadband audio source.
    let far_end_signal: Vec<f32> = (0..TOTAL_SAMPLES)
        .map(|_| rng.r#gen::<f32>() * 0.6 - 0.3) // Lower amplitude noise
        .collect();

    // Near-end signal (440Hz sine wave, from 2s to 4s) to simulate a user speaking.
    let mut near_end_signal = vec![0.0; TOTAL_SAMPLES];
    let start_sample = (SAMPLE_RATE * 2) as usize;
    let end_sample = (SAMPLE_RATE * 4) as usize;
    for i in start_sample..end_sample {
        let t = i as f32 / SAMPLE_RATE as f32;
        near_end_signal[i] = 0.5 * (2.0 * std::f32::consts::PI * 440.0 * t).sin();
    }

    // --- 2. Echo Simulation ---
    // A simple FIR filter to simulate a Room Impulse Response (RIR) with a few reflections.
    let rir = [0.6, 0.0, 0.0, -0.15, 0.0, 0.08, 0.0, 0.03];
    let mut echo_signal = vec![0.0; TOTAL_SAMPLES];
    for i in 0..TOTAL_SAMPLES {
        for (j, &coeff) in rir.iter().enumerate() {
            let delay = j * 20; // space out reflections
            if i >= delay {
                echo_signal[i] += far_end_signal[i - delay] * coeff;
            }
        }
    }

    // --- 3. Microphone Signal Creation ---
    // The final mic signal is the sum of the simulated echo and the near-end speech.
    let mic_signal: Vec<f32> = echo_signal
        .iter()
        .zip(near_end_signal.iter())
        .map(|(&echo, &near)| (echo + near).clamp(-1.0, 1.0))
        .collect();

    // --- 4. AEC Processing Loop ---
    let mut processed_signal = Vec::new();
    for i in (0..TOTAL_SAMPLES).step_by(FRAME_SIZE) {
        if i + FRAME_SIZE > mic_signal.len() || i + FRAME_SIZE > far_end_signal.len() { break; }
        
        let far_frame = &far_end_signal[i..i + FRAME_SIZE];
        let mic_frame = &mic_signal[i..i + FRAME_SIZE];

        let output_frame = aec.process(far_frame, mic_frame);
        processed_signal.extend_from_slice(&output_frame);
    }

    // --- 5. Save Results to WAV Files ---
    let output_dir = "output_generated";
    std::fs::create_dir_all(output_dir).expect("Failed to create output directory");
    
    save_wav(format!("{}/01_farend.wav", output_dir), &far_end_signal).expect("Failed to save far_end");
    save_wav(format!("{}/02_nearend.wav", output_dir), &near_end_signal).expect("Failed to save near_end");
    save_wav(format!("{}/03_mic_before_aec.wav", output_dir), &mic_signal).expect("Failed to save mic_signal");
    save_wav(format!("{}/04_mic_after_aec.wav", output_dir), &processed_signal).expect("Failed to save processed_signal");

    println!("\nSimulation complete!");
    println!("WAV files saved in '{}' directory.", output_dir);
    println!("- 03_mic_before_aec.wav: Contains echo + voice.");
    println!("- 04_mic_after_aec.wav: Should contain only voice.");
    println!("\nListen to the files to verify echo cancellation performance.");
}

/// Helper function to save a f32 signal slice to a 16-bit PCM WAV file.
fn save_wav(path: String, signal: &[f32]) -> Result<(), hound::Error> {
    let spec = hound::WavSpec {
        channels: 1,
        sample_rate: SAMPLE_RATE,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };
    let mut writer = hound::WavWriter::create(&path, spec)?;
    for &sample in signal.iter() {
        let amplitude = i16::MAX as f32;
        writer.write_sample((sample.clamp(-1.0, 1.0) * amplitude) as i16)?;
    }
    writer.finalize()?;
    Ok(())
} 