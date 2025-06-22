use nalgebra::DVector;
use num_complex::Complex;
use rustfft::{Fft, FftPlanner};
use std::sync::Arc;

/// Implements an Acoustic Echo Canceller using the Frequency Domain Adaptive Filter (FDAF)
/// algorithm with the Overlap-Save method.
///
/// This struct holds the state for the AEC and processes audio in frames.
pub struct FdafAec {
    fft_size: usize,
    frame_size: usize,
    fft: Arc<dyn Fft<f32>>,
    ifft: Arc<dyn Fft<f32>>,
    weights: DVector<Complex<f32>>,
    far_end_buffer: DVector<f32>,
    mu: f32,
    psd: DVector<f32>,
    smoothing_factor: f32,
}

impl FdafAec {
    /// Creates a new `FdafAec` instance.
    ///
    /// # Arguments
    ///
    /// * `fft_size`: The size of the FFT. This determines the filter length and the trade-off
    ///   between computational complexity and filter performance. A larger `fft_size` provides
    ///   a longer filter, which can cancel more delayed echoes, but increases latency and
    ///   computational cost. Must be a power of two.
    /// * `step_size`: The learning rate (mu) for the adaptive filter. It controls how fast the
    ///   filter adapts. A larger value leads to faster convergence but can be less stable.
    ///   A typical value is between 0.1 and 1.0.
    pub fn new(fft_size: usize, step_size: f32) -> Self {
        assert!(fft_size > 0 && fft_size.is_power_of_two(), "fft_size must be a power of two.");
        let frame_size = fft_size / 2;
        let mut fft_planner = FftPlanner::new();
        let fft = fft_planner.plan_fft_forward(fft_size);
        let ifft = fft_planner.plan_fft_inverse(fft_size);

        Self {
            fft_size,
            frame_size,
            fft,
            ifft,
            weights: DVector::from_element(fft_size, Complex::new(0.0, 0.0)),
            far_end_buffer: DVector::from_element(fft_size, 0.0),
            mu: step_size,
            psd: DVector::from_element(fft_size, 1.0), // Initialize with 1 to avoid division by zero
            smoothing_factor: 0.98,
        }
    }

    /// Processes a frame of audio data to remove echo.
    ///
    /// # Arguments
    ///
    /// * `far_end_frame`: A slice representing the audio frame from the far-end (the reference signal, e.g., loudspeaker).
    ///   Its length must be `fft_size / 2`.
    /// * `mic_frame`: A slice representing the audio frame from the near-end microphone, containing both the
    ///   near-end speaker's voice and the echo from the far-end. Its length must be `fft_size / 2`.
    ///
    /// # Returns
    ///
    /// A `Vec<f32>` containing the echo-cancelled audio frame. The length of the vector is `fft_size / 2`.
    pub fn process(&mut self, far_end_frame: &[f32], mic_frame: &[f32]) -> Vec<f32> {
        assert_eq!(far_end_frame.len(), self.frame_size, "Input far-end frame size must be half of FFT size.");
        assert_eq!(mic_frame.len(), self.frame_size, "Input mic frame size must be half of FFT size.");

        // 1. Update far-end buffer (shift old data, add new data)
        // This creates a rolling window of the last `fft_size` samples.
        self.far_end_buffer.as_mut_slice().copy_within(self.frame_size.., 0);
        self.far_end_buffer
            .rows_mut(self.frame_size, self.frame_size)
            .copy_from_slice(far_end_frame);

        // 2. FFT of the far-end signal block
        let mut x_t_buffer: Vec<Complex<f32>> = self
            .far_end_buffer
            .iter()
            .map(|&x| Complex::new(x, 0.0))
            .collect();
        self.fft.process(&mut x_t_buffer);
        let x_f = DVector::from_vec(x_t_buffer);

        // 3. Update Power Spectral Density (PSD) of the far-end signal
        for i in 0..self.fft_size {
            let power = x_f[i].norm_sqr();
            self.psd[i] = self.smoothing_factor * self.psd[i] + (1.0 - self.smoothing_factor) * power;
        }

        // 4. Estimate echo in frequency domain
        let y_f = self.weights.component_mul(&x_f);

        // 5. Inverse FFT of the estimated echo
        let mut y_t_complex = y_f.as_slice().to_vec();
        self.ifft.process(&mut y_t_complex);

        // IFFT normalization and extract real part
        let fft_size_f32 = self.fft_size as f32;
        let y_t: DVector<f32> = DVector::from_iterator(
            self.fft_size,
            y_t_complex.iter().map(|c| c.re / fft_size_f32),
        );

        // 6. Extract the valid part of the convolution (Overlap-Save method)
        let estimated_echo = y_t.rows(self.frame_size, self.frame_size);

        // 7. Calculate the error signal (mic signal - estimated echo)
        let error_signal: Vec<f32> = mic_frame
            .iter()
            .zip(estimated_echo.iter())
            .map(|(mic, echo)| mic - echo)
            .collect();

        // 8. FFT of the error signal for weight update
        // The error signal is placed in the second half of the buffer (the first half
        // is zero-padded) to ensure correct time alignment for the gradient calculation.
        let mut e_t_buffer = vec![Complex::new(0.0, 0.0); self.fft_size];
        for (i, &sample) in error_signal.iter().enumerate() {
            e_t_buffer[i + self.frame_size] = Complex::new(sample, 0.0);
        }
        
        self.fft.process(&mut e_t_buffer);
        let e_f = DVector::from_vec(e_t_buffer);
        
        // 9. Update filter weights using Normalized LMS algorithm
        let mut gradient = x_f.map(|c| c.conj()).component_mul(&e_f);
        for i in 0..self.fft_size {
            // Normalize by the PSD of the far-end signal
            gradient[i] /= self.psd[i] + 1e-10; // Add a small epsilon for stability
        }
        self.weights += &gradient * Complex::new(self.mu, 0.0);

        // 10. Return the echo-cancelled (error) signal
        error_signal
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_instance_and_process_frame() {
        const FFT_SIZE: usize = 512;
        const FRAME_SIZE: usize = FFT_SIZE / 2;
        const STEP_SIZE: f32 = 0.5;

        let mut aec = FdafAec::new(FFT_SIZE, STEP_SIZE);

        let far_end_frame = vec![0.0; FRAME_SIZE];
        let mic_frame = vec![0.1; FRAME_SIZE]; // Some non-zero value

        let error_signal = aec.process(&far_end_frame, &mic_frame);

        // Check output length
        assert_eq!(error_signal.len(), FRAME_SIZE);

        // Check for NaN or Infinity
        assert!(error_signal.iter().all(|&x| x.is_finite()), "Output contains NaN or Infinity");
    }

    #[test]
    #[should_panic]
    fn test_new_with_non_power_of_two_fft_size() {
        FdafAec::new(511, 0.5);
    }

    #[test]
    #[should_panic]
    fn test_process_with_wrong_frame_size() {
        let mut aec = FdafAec::new(512, 0.5);
        let far_end_frame = vec![0.0; 128];
        let mic_frame = vec![0.0; 256];
        aec.process(&far_end_frame, &mic_frame);
    }
}
