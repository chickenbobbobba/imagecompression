#pragma once
#include <fftw3.h>
#include <vector>
#include <complex>
#include <mutex>

class FFT {
public:
    static void init(size_t N);
    static void forward(std::vector<std::complex<double>>& data);
    static void backward(std::vector<std::complex<double>>& data);
    static void cleanup();

private:
    static fftw_plan forward_plan;
    static fftw_plan backward_plan;
    static size_t size;
    static std::mutex fftw_mutex; // FFTW is not thread-safe during plan creation
};
