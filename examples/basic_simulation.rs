//! Basic Sine Wave Simulation
//!
//! This example demonstrates the core AEC functionality in a simple, controlled environment.
//! It uses sine waves for the far-end and near-end signals and prints the Root Mean Square (RMS)
//! energy before and after processing to show the echo reduction.
//!
//! ## How to Run
//!
//! ```sh
//! cargo run --example basic_simulation
//! ```

use fdaf_aec::FdafAec;

fn main() {
    println!("--- Running Basic Sine Wave AEC Simulation ---");

    const SAMPLE_RATE: u32 = 16000;
    const FFT_SIZE: usize = 1024;
    const FRAME_SIZE: usize = FFT_SIZE / 2;
    const STEP_SIZE: f32 = 0.05;

    let mut aec = FdafAec::new(FFT_SIZE, STEP_SIZE);

    // --- Signal Generation ---
    // 1. Far-end signal: A 440Hz sine wave, representing the audio from the loudspeaker.
    let mut far_end_signal = Vec::new();
    for i in 0..(SAMPLE_RATE * 2) { // 2 seconds of audio
        let t = i as f32 / SAMPLE_RATE as f32;
        far_end_signal.push(0.6 * (2.0 * std::f32::consts::PI * 440.0 * t).sin());
    }

    // 2. Near-end signal: An 880Hz sine wave, representing the local user's voice.
    //    It starts after 0.5 seconds to create a period of single-talk (echo only).
    let mut near_end_signal = vec![0.0; (SAMPLE_RATE * 2) as usize];
    for i in (SAMPLE_RATE / 2) as usize..near_end_signal.len() {
        let t = i as f32 / SAMPLE_RATE as f32;
        near_end_signal[i] = 0.4 * (2.0 * std::f32::consts::PI * 880.0 * t).sin();
    }

    // 3. Microphone signal: A mix of the near-end signal and a delayed, attenuated
    //    version of the far-end signal (the echo).
    let echo_delay_samples = 128;
    let echo_attenuation = 0.7;
    let mut mic_signal = vec![0.0; far_end_signal.len()];
    for i in echo_delay_samples..mic_signal.len() {
        mic_signal[i] = near_end_signal[i] + far_end_signal[i - echo_delay_samples] * echo_attenuation;
    }

    // --- AEC Processing Loop ---
    let mut processed_signal = Vec::new();
    for (mic_chunk, far_chunk) in mic_signal.chunks(FRAME_SIZE).zip(far_end_signal.chunks(FRAME_SIZE)) {
        if mic_chunk.len() != FRAME_SIZE { break; }

        let output_frame = aec.process(far_chunk, mic_chunk);
        processed_signal.extend_from_slice(&output_frame);
    }

    // --- Performance Analysis ---
    println!("\n--- AEC Performance Analysis (RMS Energy) ---");
    analyze_and_print_rms("Single-Talk (Echo Only)",
        &mic_signal[.. (SAMPLE_RATE / 2) as usize],
        &processed_signal[.. (SAMPLE_RATE / 2) as usize]);

    analyze_and_print_rms("Double-Talk (Echo + Voice)",
        &mic_signal[(SAMPLE_RATE / 2) as usize ..],
        &processed_signal[(SAMPLE_RATE / 2) as usize ..]);

    println!("\nExplanation:");
    println!("- Single-Talk: The RMS of the processed signal should be significantly lower, showing echo removal.");
    println!("- Double-Talk: The RMS should decrease (as echo is removed) but not go to zero, showing the near-end voice was preserved.");
}

/// Helper function to calculate RMS for a signal slice.
fn rms(signal: &[f32]) -> f32 {
    if signal.is_empty() { return 0.0; }
    let sum_sq: f32 = signal.iter().map(|&x| x * x).sum();
    (sum_sq / signal.len() as f32).sqrt()
}

/// Helper function to print analysis for a signal segment.
fn analyze_and_print_rms(label: &str, before: &[f32], after: &[f32]) {
    println!("\n[{}]", label);
    println!(" - Before AEC: {:.6}", rms(before));
    println!(" - After AEC:  {:.6}", rms(after));
} 