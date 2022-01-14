use core::f64::consts::PI;
use num::Complex;
use num_traits::cast;
use num_traits::float::{Float, FloatConst};
use num_traits::identities::{one, zero};

pub struct Fft<T> {
    // フーリエ変換結果の次元数
    n: usize,
    // ビットリバースしたインデックスを付与
    r_indexes: Vec<usize>,
    // 事前計算したW_N^kを保持
    w: Vec<Complex<T>>,
}

// Float + FloatConst でfloat型を許可
impl<T: Float + FloatConst + std::fmt::Debug> Fft<T> {
    pub fn new() -> Self {
        Self {
            n: 0,
            r_indexes: vec![],
            w: vec![],
        }
    }

    pub fn setup(&mut self, n: usize) {
        let n_power2 = n.trailing_zeros();
        // nは2の累乗である前提
        assert_eq!(n, 1 << n_power2, "len of n should be 2^x");

        self.n = n;
        self.r_indexes = self.calc_bit_reversed_indexes(n);
        self.w = self.calc_w(n);
    }

    fn calc_bit_reversed_indexes(&mut self, n: usize) -> Vec<usize> {
        let n_power2 = n.trailing_zeros();
        // このメソッドは、nが2の累乗である前提で作成
        assert_eq!(n, 1 << n_power2, "len of n should be 2^x");

        let mut r_indexes = Vec::<usize>::with_capacity(n);

        // 最初に /0b0*10*/ を代入 ((power+1)番目のビット数が1)
        let mut r_bit = 1 << n_power2;

        // n_power2番目以下のビットを扱う
        // 2bitを1単位として扱い、最後に2が残る場合1bitで演算

        // e.g. n_power2 = 5, r_bit = 2^5 = 32
        //  [r_bit = 32 = 0b100000]
        //    r[ 0] = r[0b00000] = 0b00000 ++
        //  [r_bit = 8 = 0b001000]
        //    r[ 0] = r[0b00000] = 0b00000
        //    r[ 1] = r[0b000(01)] = 0b(01)000 ++
        //    r[ 2] = r[0b000(10)] = 0b(10)000 ++
        //    r[ 3] = r[0b000(11)] = 0b(11)000 ++
        //  [r_bit = 2 = 0b000(01)0]
        //    r[ 0] = r[0b00000] = 0b00000
        //    r[ 1] = r[0b000(01)] = 0b(01)000
        //    r[ 2] = r[0b000(10)] = 0b(10)000
        //    r[ 3] = r[0b000(11)] = 0b(11)000
        //    r[ 4] = r[0b0(01)(00)] = 0b(00)(01)0 ++
        //    r[ 5] = r[0b0(01)(01)] = 0b(01)(01)0 ++
        //    r[ 6] = r[0b0(01)(10)] = 0b(10)(01)0 ++
        //    r[ 7] = r[0b0(01)(11)] = 0b(11)(01)0 ++
        //    r[ 8] = r[0b0(10)(00)] = 0b(00)(10)0 ++
        //    r[ 9] = r[0b0(10)(01)] = 0b(01)(10)0 ++
        //     :
        //    r[15] = r[0b0(11)(11)] = 0b(11)(11)0 ++
        //  [r_bit = 1 = 0b000001]
        //    r[ 0] = r[0b00000] = 0b00000
        //    r[ 1] = r[0b000(01)] = 0b(01)000
        //     :
        //    r[15] = r[0b0(11)(11)] = 0b(11)(11)0
        //    r[16] = r[0b(1)(00)(00)] = 0b(00)(00)(1) ++
        //    r[17] = r[0b(1)(00)(01)] = 0b(01)(00)(1) ++
        //     :
        //    r[31] = r[0b(1)(11)(11)] = 0b(11)(11)(1) ++

        // 最初に /0b0*/ を追加
        r_indexes.push(0);
        // 4の累乗部として計算
        while r_bit > 2 {
            // ビット数を2つ右にずらす
            r_bit >>= 2;
            // 2bitを1単位とし、(01),(10),(11)を追加
            let len = r_indexes.len();
            for j in 0..len {
                r_indexes.push(r_indexes[j] | r_bit);
            }
            for j in 0..len {
                r_indexes.push(r_indexes[j] | r_bit << 1);
            }
            for j in 0..len {
                r_indexes.push(r_indexes[j] | r_bit | r_bit << 1);
            }
        }
        // 2が残る場合、最後に基底数2でバタフライ演算を行うため、indexも最後に算出
        if r_bit == 2 {
            for j in 0..r_indexes.len() {
                r_indexes.push(r_indexes[j] | 1);
            }
        }

        // in-place演算用にインデックスを加工
        return self.convert_indexes_as_inplace(r_indexes);
    }

    fn convert_indexes_as_inplace(&mut self, r_indexes: Vec<usize>) -> Vec<usize> {
        let mut nums = (0..r_indexes.len()).collect::<Vec<_>>();

        // 整数列のインデックスへr_indexesでswapを実行
        // その際のインデックスを保持して実データのswapに利用することで、結果的にr_indexesで整頓したことになる
        return (0..r_indexes.len())
            .map(|i| {
                let r_i = r_indexes[i];
                let swapped_r_i = (0..nums.len()).find(|&j| nums[j] == r_i).unwrap();
                nums.swap(i, swapped_r_i);
                swapped_r_i
            })
            .collect();
    }

    fn calc_w(&mut self, n: usize) -> Vec<Complex<T>> {
        let n_power2 = n.trailing_zeros();
        // このメソッドは、nが2の累乗である前提で作成
        assert_eq!(n, 1 << n_power2, "len of n should be 2^x");

        let mut w = Vec::with_capacity(n + 1);

        // W_N^0=1
        w.push(one());

        // nが2以下の場合
        if n <= 2 {
            if n == 2 {
                // W_N^N/2=-1
                w.push(cast(-1.0).unwrap());
            }
            // W_N^N=1
            w.push(one());
            return w;
        }

        // nが2以下ではない場合(=nが4で割り切れる場合)
        let q = n >> 2;
        let h = n >> 1;
        // 0~N/4(実直に計算)
        for i in 1..q {
            w.push(self.calc_part_w(n, i));
        }
        // W_N^N/4=-i
        w.push(-Complex::i());
        // N/4～N/2(計算結果を流用)
        for i in q + 1..h {
            let tmp = w[i - q];
            w.push(Complex::new(tmp.im, -tmp.re));
        }
        // W_N^N/2=-1
        w.push(cast(-1.0).unwrap());
        // N/2～N(計算結果を流用)
        for i in h + 1..n {
            let tmp = w[i - h];
            w.push(Complex::new(-tmp.re, -tmp.im));
        }
        // W_N^N=1
        w.push(one());

        return w;
    }

    fn calc_part_w(&mut self, n: usize, seq: usize) -> Complex<T> {
        // e^(-i2π)をN分割。seqごとの値を取得
        Complex::new(
            zero(),
            cast::<_, T>(-2.0).unwrap() * T::PI() / cast(n).unwrap() * cast(seq).unwrap(),
        )
        .exp()
    }

    pub fn process(&mut self, frames: &mut [Complex<T>]) {
        let len = frames.len();
        let len_power2 = len.trailing_zeros();
        // このメソッドは、フレーム数が2の累乗である前提で作成
        assert_eq!(len, 1 << len_power2, "len of frames should be 2^x");

        // 1フレーム以下の場合、そのまま返す
        if len <= 1 {
            return;
        }

        // フレーム数がnと異なる場合、nをlenに置き換える
        if len != self.n {
            self.setup(len);
        }

        self.inner_process(frames);
    }

    // 処理結果を格納するVectorはin-palce演算なので直接もらう
    fn inner_process(&mut self, x: &mut [Complex<T>]) {
        let n = x.len();

        // whileループごとのインデックス数(x_kにおけるkを決める)
        let mut calc_n = n;
        let mut calc_w_bit = 0;
        // 4 or 2になるまでバタフライ演算(2bit)
        while calc_n > 4 {
            let before_n = calc_n;
            calc_n >>= 2;
            for i in 0..calc_n {
                let w_i = i << calc_w_bit;
                let (w1, w2, w3) = (self.w[w_i], self.w[w_i << 1], self.w[w_i * 3]);
                for begin_i in (0..n).step_by(before_n) {
                    // 各インデックスにはcalc_n分の間が空く
                    let i0 = begin_i + i;
                    let i1 = i0 + calc_n;
                    let i2 = i1 + calc_n;
                    let i3 = i2 + calc_n;
                    // xの値を計算
                    // 以下の計算を短縮し高速化
                    // let x0 = x[i0] + x[i1] + x[i2] + x[i3];
                    // let x1 = (x[i0] = i() * x[i1] - x[i2] + i() * x[i3]) * w1;
                    // let x2 = (x[i0] - x[i1] + x[i2] - x[i3]) * w2;
                    // let x3 = (x[i0] + i() * x[i1] - x[i2] - i() * x[i3]) * w3;
                    let xi0_plus_xi2 = x[i0] + x[i2];
                    let xi0_minus_xi2 = x[i0] - x[i2];
                    let xi1_plus_xi3 = x[i1] + x[i3];
                    let xi1_minus_xi3 = x[i1] - x[i3];
                    let xi1_minus_xi3_i = Complex::new(-xi1_minus_xi3.im, xi1_minus_xi3.re);
                    // 変形前の値を上書き(=in-place演算)
                    x[i0] = xi0_plus_xi2 + xi1_plus_xi3;
                    x[i1] = (xi0_minus_xi2 - xi1_minus_xi3_i) * w1;
                    x[i2] = (xi0_plus_xi2 - xi1_plus_xi3) * w2;
                    x[i3] = (xi0_minus_xi2 + xi1_minus_xi3_i) * w3;
                }
            }
            calc_w_bit += 2;
        }

        // 最後のバタフライ演算は処理を簡略化
        if calc_n == 4 {
            for i0 in (0..n).step_by(calc_n) {
                let i1 = i0 + 1;
                let i2 = i1 + 1;
                let i3 = i2 + 1;
                // xの値を計算
                // 以下の計算を短縮し高速化
                // let x0 = x[i0] + x[i1] + x[i2] + x[i3];
                // let x1 = (x[i0] = i() * x[i1] - x[i2] + i() * x[i3]) * w1;
                // let x2 = (x[i0] - x[i1] + x[i2] - x[i3]) * w2;
                // let x3 = (x[i0] + i() * x[i1] - x[i2] - i() * x[i3]) * w3;
                let xi0_plus_xi2 = x[i0] + x[i2];
                let xi0_minus_xi2 = x[i0] - x[i2];
                let xi1_plus_xi3 = x[i1] + x[i3];
                let xi1_minus_xi3 = x[i1] - x[i3];
                let xi1_minus_xi3_i = Complex::new(-xi1_minus_xi3.im, xi1_minus_xi3.re);
                // 変形前の値を上書き(=in-place演算)
                x[i0] = xi0_plus_xi2 + xi1_plus_xi3;
                x[i1] = xi0_minus_xi2 - xi1_minus_xi3_i;
                x[i2] = xi0_plus_xi2 - xi1_plus_xi3;
                x[i3] = xi0_minus_xi2 + xi1_minus_xi3_i;
            }
        } else {
            // calc_n = 2
            for i in (0..n).step_by(calc_n) {
                let even_i = i;
                let odd_i = i + 1;
                // xの値を計算
                let x_even = x[even_i] + x[odd_i];
                let x_odd = x[even_i] - x[odd_i];
                // 変形前の値を上書き(=in-place演算)
                x[even_i] = x_even;
                x[odd_i] = x_odd;
            }
        }

        // in-place演算として高速化するため、ビット反転をメモリスワップで行う
        for (i, &r) in self.r_indexes.iter().enumerate() {
            x.swap(i, r);
        }
    }
}

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
    assert_eq!(n, 1 << n_power2);

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
    let mut r_bit = 1 << power;

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
