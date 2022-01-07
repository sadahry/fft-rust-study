use fft_rs::dtft;
use std::f64::consts::PI;

fn main() {
    // create sin_curve fn
    let flame_len = 1000;
    let samplerate = 16000;
    let hz = 2000;
    let sin_curve: Vec<f64> = (0..flame_len)
        .map(|a| (a as f64) * 2.0 * PI * (hz as f64) / (samplerate as f64))
        .map(|a| a.sin())
        .collect();

    // ./target/release/fft_rs  2.85s user 0.02s system 97% cpu 2.938 total
    let rslt = dtft(sin_curve, samplerate);

    let (max_index, max) =
        rslt.iter()
            .enumerate()
            .fold((usize::MIN, f64::MIN), |(i_a, a), (i_b, &b)| {
                if b.norm() > a {
                    (i_b, b.norm())
                } else {
                    (i_a, a)
                }
            });
    println!("rslt len: {:?}", rslt.len()); // = fs/2
    println!("max rslt f: {:?}", max_index + 1); // 1 <= f <= fs/2
    println!("max rslt: {:?}", max);
}
