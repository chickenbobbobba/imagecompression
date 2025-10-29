#include <vector>
#include <complex>
#include <cassert>

bool is_pow_two(int a) {
    while (a > 1) {
        if (a % 2 == 1) return false; 
        a /= 2;
    }
    return true;
}

long revbits(const long a, const long befbit) {
    long b = 0; 
    for (int i = 0; i < befbit; i++) {
        b <<= 1;
        b += (a & (1<<i)) != 0;
    }
    return b;
}

std::vector<std::complex<double>> FFT_ip(const std::vector<std::complex<double>>& data) {
    assert(is_pow_two(data.size()) == true);

    size_t n = data.size();
    int rotpoint = std::countr_zero(n);
    
    std::vector<std::complex<double>> result = data;

    for (size_t i = 0; i < n; i++) {
        result[revbits(i, rotpoint)] = data[i];
    }
    
    std::vector<std::complex<double>> twiddles(n / 2);
    for (size_t k = 0; k < n / 2; ++k)
        twiddles[k] = std::polar(1.0, -2.0 * M_PI * k / n);

    for (int len = 2; len <= n; len *= 2) {
        for (int batchPtr = 0; batchPtr < n; batchPtr += len) {
            int half = len/2;
            size_t step = n / len;

            for (int j = 0; j < half; j++) {
                int evenIdx = batchPtr + j;
                int oddIdx  = evenIdx + half;

                auto twiddle = twiddles[j * step];

                auto u = result[evenIdx];
                auto v = result[oddIdx] * twiddle;

                result[evenIdx] = u + v;
                result[oddIdx]  = u - v;
            }
        }
    }
    
    return result;
}

std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>& data);

std::vector<std::complex<double>> IFFT(const std::vector<std::complex<double>>& data) {
    std::vector<std::complex<double>> conjd(data.size());
    for (size_t i = 0; i < data.size(); ++i) conjd[i] = std::conj(data[i]);

    auto res = FFT(conjd);
    for (auto& v : res) v = std::conj(v) / static_cast<double>(data.size());
    return res;
}

std::vector<std::complex<double>> IFFT_ip(const std::vector<std::complex<double>>& data) {
    std::vector<std::complex<double>> conjd(data.size());
    for (size_t i = 0; i < data.size(); ++i) conjd[i] = std::conj(data[i]);

    auto res = FFT_ip(conjd);
    for (auto& v : res) v = std::conj(v) / static_cast<double>(data.size());
    return res;
}

std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>& data) {
    const size_t n = data.size();
    if (n == 0) return {};

    if (is_pow_two(n)) return FFT_ip(data);

    size_t m = 1;
    while (m < 2 * n - 1) m <<= 1;

    const double PI = std::acos(-1.0);

    std::vector<std::complex<double>> a(m, 0.0), b(m, 0.0);

    for (size_t i = 0; i < n; ++i) {
        double angle = PI * (i * i) / (double)n;
        std::complex<double> w = std::polar(1.0, -angle);
        a[i] = data[i] * std::conj(w);
        b[i] = w;
        if (i != 0) b[m - i] = w;
    }

    auto A = FFT_ip(a);
    auto B = FFT_ip(b);
    for (size_t i = 0; i < m; ++i)
        A[i] *= B[i];
    A = IFFT_ip(A);

    std::vector<std::complex<double>> result(n);
    for (size_t i = 0; i < n; ++i) {
        double angle = PI * (i * i) / n;
        std::complex<double> w(std::cos(angle), std::sin(angle));
        result[i] = A[i] * w;
    }

    return result;
}