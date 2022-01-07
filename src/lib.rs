use num::complex::Complex;
use std::f64::consts::PI;

pub fn dtft(frames: Vec<f64>, fs: i64) -> Vec<Complex<f64>> {
    let mut rslt: Vec<Complex<f64>> = Vec::new();
    // 1 <= f <= fs/2
    for f in 1..=(fs / 2) {
        let mut sigma: Complex<f64> = Complex::new(0.0, 0.0);
        // 0 <= k <= N-1
        for (k, xk) in frames.iter().enumerate() {
            // t = k / fs
            let t: f64 = (k as f64) / (fs as f64);
            // ω = 2πf
            let w = 2.0 * PI * (f as f64);
            // e^(-iωt)
            let exp = Complex::new(0.0, -w * t).exp();
            // X(f) = ∑xk*e^(-iωt)
            sigma += xk * exp;
        }
        rslt.push(sigma);
    }
    return rslt;
}
