#include <fftw3.h>
#include <fftwrap.hpp>
#include <cassert>

fftw_plan FFT::forward_plan = nullptr;
fftw_plan FFT::backward_plan = nullptr;
size_t FFT::size = 0;
std::mutex FFT::fftw_mutex;

void FFT::init(size_t N) {
    std::lock_guard<std::mutex> lock(fftw_mutex);
    if (N == size && forward_plan && backward_plan) return;

    cleanup(); // free any old plans
    size = N;

    std::vector<std::complex<double>> temp(N);

    forward_plan = fftw_plan_dft_1d(
        N,
        reinterpret_cast<fftw_complex*>(temp.data()),
        reinterpret_cast<fftw_complex*>(temp.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );

    backward_plan = fftw_plan_dft_1d(
        N,
        reinterpret_cast<fftw_complex*>(temp.data()),
        reinterpret_cast<fftw_complex*>(temp.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    fftw_print_plan(forward_plan);
}

void FFT::forward(std::vector<std::complex<double>>& data) {
    assert(data.size() == size);
    fftw_execute_dft(forward_plan,
        reinterpret_cast<fftw_complex*>(data.data()),
        reinterpret_cast<fftw_complex*>(data.data()));
}

void FFT::backward(std::vector<std::complex<double>>& data) {
    assert(data.size() == size);
    fftw_execute_dft(backward_plan,
        reinterpret_cast<fftw_complex*>(data.data()),
        reinterpret_cast<fftw_complex*>(data.data()));
}

void FFT::cleanup() {
    if (forward_plan) fftw_destroy_plan(forward_plan);
    if (backward_plan) fftw_destroy_plan(backward_plan);
    forward_plan = backward_plan = nullptr;
}