#pragma once
#include <fftw3.h>
#include <vector>
#include <complex>
#include <mutex>
#include <unordered_map>

class FFT {
private:
    struct PlanPair;  // Forward declaration
    static std::unordered_map<size_t, PlanPair> plans;
    static std::mutex fftw_mutex;
    
public:
    static void init(size_t N);
    static void forward(std::vector<std::complex<double>>& data);
    static void backward(std::vector<std::complex<double>>& data);
    static void cleanup();
};