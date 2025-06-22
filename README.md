# FdafAec: Frequency Domain Adaptive Filter for Acoustic Echo Cancellation

[![Crates.io](https://img.shields.io/crates/v/fdaf-aec.svg)](https://crates.io/crates/fdaf-aec)
[![Docs.rs](https://docs.rs/fdaf-aec/badge.svg)](https://docs.rs/fdaf-aec)

`fdaf-aec` is a Rust library that implements an Acoustic Echo Canceller (AEC) using the Frequency Domain Adaptive Filter (FDAF) algorithm with the Overlap-Save method. It is designed to be simple to use while providing a robust solution for real-time echo cancellation in voice communication applications.

## What is Acoustic Echo Cancellation?

In a typical hands-free communication setup (like a video call), the sound from the far-end user (played through your speakers) can be picked up by your microphone. This causes the far-end user to hear an echo of their own voice, which is disruptive. An AEC is a signal processing component that identifies and removes this unwanted echo from the microphone signal, resulting in a clean audio stream.

![AEC Diagram](https://i.imgur.com/GCRs5pE.png)

## How It Works

This library uses a Frequency Domain Adaptive Filter (FDAF). Here's a simplified overview of the process:

1.  **Buffering**: It takes a frame of audio from the far-end (speaker) signal and the microphone signal.
2.  **FFT**: It transforms these audio signals into the frequency domain using the Fast Fourier Transform (FFT).
3.  **Echo Estimation**: In the frequency domain, an adaptive filter (represented by a set of complex weights) models the echo path. It uses the far-end signal to predict what the echo should sound like.
4.  **Subtraction**: The predicted echo is subtracted from the microphone signal, leaving (ideally) only the near-end user's voice.
5.  **Adaptation**: The filter constantly adjusts its weights using the Normalized Least Mean Squares (NLMS) algorithm to adapt to changing room acoustics and echo paths.
6.  **IFFT**: The cleaned signal is transformed back into the time domain (audio samples) and returned.

The **Overlap-Save** method is used to efficiently process the audio in blocks, making it suitable for real-time applications.

## Features

- Real-time capable FDAF implementation.
- Adjustable learning rate (step size) to balance convergence speed and stability.
- Simple and straightforward API.
- Minimal dependencies for the core library.

## Getting Started

To use `fdaf-aec` in your project, add it to your `Cargo.toml`:

```toml
[dependencies]
fdaf-aec = "0.1.0" # Replace with the latest version
```

Then, initialize the `FdafAec` and process your audio frames:

```rust
use fdaf_aec::FdafAec;

// Configuration
const FFT_SIZE: usize = 1024; // Must be a power of two. Determines filter length.
const FRAME_SIZE: usize = FFT_SIZE / 2; // Should be half of FFT_SIZE.
const STEP_SIZE: f32 = 0.02; // Learning rate. A small value between 0 and 1.

// Create a new AEC instance
let mut aec = FdafAec::new(FFT_SIZE, STEP_SIZE);

// Your audio processing loop
// loop {
    // let far_end_frame: &[f32] = get_far_end_audio_frame(); // Should have FRAME_SIZE samples
    // let mic_frame: &[f32] = get_mic_audio_frame();     // Should have FRAME_SIZE samples

    // Process the frames to get the echo-cancelled signal
    // let output_frame = aec.process(far_end_frame, mic_frame);

    // Use the `output_frame` for your application...
// }
```

## Examples

The project includes several examples in the `examples/` directory to demonstrate its functionality.

### 1. Basic Simulation

This example runs a simple simulation with sine waves to show the core principle of echo reduction. It prints the RMS energy before and after processing.

```sh
cargo run --example basic_simulation
```

### 2. Generated Signal AEC

This example provides a more realistic test by generating white noise (far-end) and a simulated voice signal (near-end). It models a simple room echo and saves the audio before and after AEC to WAV files in the `output_generated/` directory for listening.

```sh
# The --release flag is recommended for faster execution
cargo run --example generated_signal_aec --release
```

### 3. File-Based AEC (CLI Tool)

This is a command-line utility for processing your own WAV files. It's the most practical way to test the AEC's performance on real-world recordings.

**Usage:**

```sh
# Place your WAV files in the project root directory
# For best results, use mono, 16kHz WAV files.
cargo run --example file_based_aec --release -- \
  --farend your_farend_file.wav \
  --mic your_mic_file.wav \
  --output processed_output.wav
```

## License

This project is licensed under the MIT License.
