use fft_rs::*;
use num::Complex;
use std::f64::consts::PI;
use std::time::Instant;

fn main() {
    // create sin_curve fn
    let flame_len: usize = 512;
    let samplerate = 16000;
    let hz = 2400;
    let mut sin_curve: Vec<Complex<f64>> = (0..flame_len)
        .map(|a| (a as f64) * 2.0 * PI * (hz as f64) / (samplerate as f64))
        .map(|a| Complex::new(a.sin(), 0.0))
        .collect();

    let size = flame_len.next_power_of_two();
    let mut fft = Fft::<f64>::new();
    fft.setup(size);

    // calc time
    let start = Instant::now();
    // let rslt = dtft(&sin_curve, size);
    // let rslt = cooley_tukey_fft(&sin_curve);
    fft.process(&mut sin_curve);
    let rslt = sin_curve;
    let end = start.elapsed();

    let time = (end.as_secs() * 1000_000_000) + (end.subsec_nanos() as u64);
    println!("time: {:?}", time);

    let spectrum = &rslt[..(size / 2 as usize)];

    let (max_index, max) =
        spectrum
            .iter()
            .enumerate()
            .fold((usize::MIN, f64::MIN), |(i_a, a), (i_b, &b)| {
                if b.norm() > a {
                    (i_b, b.norm())
                } else {
                    (i_a, a)
                }
            });
    let max_hz: f64 = (((max_index) * samplerate) as f64) / (size as f64);

    println!("flame_len: {:?}", flame_len); // = 400;
    println!("rslt len: {:?}", rslt.len()); // = 256 (= 512/2);
    println!("max index: {:?}", max_index);
    println!("max rslt: {:?}", max);
    // 高速化のため、正確なヘルツとはならない
    // (samplerate / size)Hz ほどの誤差が発生する
    println!("max hz: {:?}", max_hz); // = 2000 - (31+1) <= max_hz <= 2000 + (31+1)
}
