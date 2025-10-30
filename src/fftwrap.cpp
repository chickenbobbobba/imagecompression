#include <fftw3.h>
#include <fftwrap.hpp>
#include <cassert>
#include <iostream>
#include <unordered_map>

struct FFT::PlanPair {
    fftw_plan forward;
    fftw_plan backward;
};

std::unordered_map<size_t, FFT::PlanPair> FFT::plans;
std::mutex FFT::fftw_mutex;

void FFT::init(size_t N) {
    std::lock_guard<std::mutex> lock(fftw_mutex);
    
    // Check if plan already exists for this size
    if (plans.find(N) != plans.end()) {
        return;
    }
    
    // Create new plans for this size
    std::vector<std::complex<double>> temp(N);
    
    fftw_plan forward = fftw_plan_dft_1d(
        N,
        reinterpret_cast<fftw_complex*>(temp.data()),
        reinterpret_cast<fftw_complex*>(temp.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );
    
    fftw_plan backward = fftw_plan_dft_1d(
        N,
        reinterpret_cast<fftw_complex*>(temp.data()),
        reinterpret_cast<fftw_complex*>(temp.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );
    
    plans[N] = {forward, backward};
    
    fftw_print_plan(forward);
    std::cout << "\n";
}

void FFT::forward(std::vector<std::complex<double>>& data) {
    size_t N = data.size();
    
    auto it = plans.find(N);
    assert(it != plans.end() && "FFT::init must be called before forward");
    
    fftw_execute_dft(it->second.forward,
        reinterpret_cast<fftw_complex*>(data.data()),
        reinterpret_cast<fftw_complex*>(data.data()));
}

void FFT::backward(std::vector<std::complex<double>>& data) {
    size_t N = data.size();
    
    auto it = plans.find(N);
    assert(it != plans.end() && "FFT::init must be called before backward");
    
    fftw_execute_dft(it->second.backward,
        reinterpret_cast<fftw_complex*>(data.data()),
        reinterpret_cast<fftw_complex*>(data.data()));
}

void FFT::cleanup() {
    std::lock_guard<std::mutex> lock(fftw_mutex);
    
    for (auto& [size, plan_pair] : plans) {
        if (plan_pair.forward) fftw_destroy_plan(plan_pair.forward);
        if (plan_pair.backward) fftw_destroy_plan(plan_pair.backward);
    }
    
    plans.clear();
}