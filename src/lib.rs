use num::Complex;
use std::f64::consts::PI;

pub fn dtft(frames: &[Complex<f64>], size: usize) -> Vec<Complex<f64>> {
    let mut rslt: Vec<Complex<f64>> = Vec::new();
    // 1 <= f <= size
    for f in 1..=(size) {
        let mut sigma: Complex<f64> = Complex::new(0.0, 0.0);
        // 0 <= k <= N-1
        for (k, xk) in frames.iter().enumerate() {
            // t = k / size
            let t: f64 = (k as f64) / (size as f64);
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

pub fn cooley_tukey_fft(frames: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let n = frames.len();
    // 先頭bitから2進数中のゼロを取得
    // nが2の累乗の場合"ゼロの個数＝2の累乗時の指数"となる
    let n_power2 = n.trailing_zeros();
    // nが2の累乗であることを確認
    assert_eq!(n, 2usize.pow(n_power2));

    // 処理結果を格納するVector
    // 最初は入力フレームの値(x_k)を代入
    let mut x = frames.to_vec();

    // whileループごとのインデックス数(x_kにおけるkを決める)
    let mut calc_n = n;
    // Wにおける指数の固定値を保持
    let mut w_power = -2.0 * PI / (n as f64);
    while calc_n > 1 {
        let before_n = calc_n;
        calc_n >>= 1;
        for i in 0..calc_n {
            // W_N^kの算出
            let w = Complex::new(0.0, w_power * (i as f64)).exp();
            // バタフライ演算(=偶数部と奇数部の演算)
            //// whileループごとに処理パターンが増える(変形ごとに集約されたxを利用するため)
            //// e.g. n = 8
            ////  [calc_n = 4]
            ////    1パターン: x_k
            ////    i = k, 0 <= k <= N/2-1
            ////  [calc_n = 2]
            ////    2パターン: x^(2m), x^(2m+1)
            ////    i = m, 0 <= m <= N/4-1
            ////  [calc_n = 1]
            ////    4パターン: x^(2(2l)), x^(2(2l+1)+1), x^(2(2l)), x^(2(2l+1)+1)
            ////    i = l, 0 <= l <= N/8-1
            for begin_i in (0..n).step_by(before_n) {
                // 偶数部と奇数部のインデックスにはcalc_n分の間が空く
                let even_i = begin_i + i;
                let odd_i = begin_i + i + calc_n;
                // xの値を計算
                let x_even = x[even_i] + x[odd_i];
                let x_odd = (x[even_i] - x[odd_i]) * w;
                // 変形前の値を上書き(=in-place演算)
                x[even_i] = x_even;
                x[odd_i] = x_odd;
            }
        }
        w_power *= 2.0;
    }

    // xのindexをlog2(n)ビットでビット反転することで正しい順序を取得
    return bit_reversed_indexes(n_power2)
        .iter()
        .map(|&r_i| x[r_i])
        .collect();
}

fn bit_reversed_indexes(power: u32) -> Vec<usize> {
    let mut r_indexes: Vec<usize> = Vec::new();

    // 最初に /0b0*10*/ を代入 ((power+1)番目のビット数が1)
    let mut r_bit = 2usize.pow(power);

    // power番目以下のビットを扱う
    // 最初に /0b0*/ を代入
    r_indexes.push(0);
    while r_bit > 1 {
        // ビット数を1つ右にずらす
        r_bit >>= 1;
        for j in 0..r_indexes.len() {
            // whileループごとに、処理中の配列にr_bitを加算した値をすべてpushすることで、
            // indexに対してpower番目までのビットを反転させたindex配列が作成される
            // e.g. power = 3, r_bit = 2^3 = 8
            //  [r_bit = 8 = 0b1000]
            //    r[0] = r[0b000] = 0b000 ++
            //  [r_bit = 4 = 0b0100]
            //    r[0] = r[0b000] = 0b000
            //    r[1] = r[0b001] = 0b100 ++
            //  [r_bit = 2 = 0b0010]
            //    r[0] = r[0b000] = 0b000
            //    r[1] = r[0b001] = 0b100
            //    r[2] = r[0b010] = 0b010 ++
            //    r[3] = r[0b011] = 0b110 ++
            //  [r_bit = 1 = 0b0001]
            //    r[0] = r[0b000] = 0b000
            //    r[1] = r[0b001] = 0b100
            //    r[2] = r[0b010] = 0b010
            //    r[3] = r[0b011] = 0b110
            //    r[0] = r[0b100] = 0b001 ++
            //    r[1] = r[0b101] = 0b101 ++
            //    r[2] = r[0b110] = 0b011 ++
            //    r[3] = r[0b111] = 0b111 ++
            r_indexes.push(r_indexes[j] | r_bit);
        }
    }
    return r_indexes;
}
