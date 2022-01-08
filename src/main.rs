use fft_rs::dtft;
use std::f64::consts::PI;

fn get_size(flame_len: i64) -> i64 {
    let mut i: u32 = 0;
    while flame_len > 2i64.pow(i) {
        i += 1;
    }
    2i64.pow(i)
}

fn main() {
    // create sin_curve fn
    let flame_len = 400;
    let samplerate = 16000;
    let hz = 2000;
    let sin_curve: Vec<f64> = (0..flame_len)
        .map(|a| (a as f64) * 2.0 * PI * (hz as f64) / (samplerate as f64))
        .map(|a| a.sin())
        .collect();

    let size = get_size(flame_len);
    let rslt = dtft(sin_curve, size);

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
    let max_hz: f64 = (((max_index + 1) * samplerate) as f64) / (size as f64);

    println!("flame_len: {:?}", flame_len); // = 400;
    println!("rslt len: {:?}", rslt.len()); // = 256 (= 512/2);
    println!("max index: {:?}", max_index);
    println!("max rslt: {:?}", max);
    println!("max hz: {:?}", max_hz); // = 2000
}
