#include <algorithm>
#include <complex>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <hilbert.h>
#include <fstream>
#include <fft.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

void savePPM(const std::string& path, int width, int height, const std::vector<unsigned char>& data) {
    std::cout << width * height * 3 << " " << data.size() << "\n";
    std::ofstream out(path, std::ios::binary);
    out << "P6\n" << width << " " << height << "\n255\n";
    out.write((char*)data.data(), data.size());
}

int main(int argc, char** argv){
    double keepfrac = /* 1/ */1;
    for (;;) {
        std::string filepath = "../resources/ball.png";
        if (argc > 1) filepath = argv[1];
        std::cout << "filepath: " << filepath << "\n";
        std::cout << "reading...\n";
        int width = 0;
        int height = 0;
        int channels = 0;
        
        unsigned char* rawIn = stbi_load(filepath.c_str(), &width, &height, &channels, 3);
        std::cout << "width: " << width << "\nheight: " << height << "\nchannels: " << channels << "\n";

        std::cout << "encoding...\n";

        std::vector<unsigned char> rawOut(3 * width * height);
        
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int index = 3 * (j + width * i);
                int hilbidx = 3 * gilbidx(j, i, width, height);

                for (int k = 0; k < channels; k++) {
                    rawOut[hilbidx + k] = rawIn[index + k];
                }
            }
        }

        std::cout << "FFT prep...\n";
        std::vector<std::complex<double>> waves;
        waves.reserve(rawOut.size());

        for (auto i : rawOut) {
            waves.emplace_back(i, 0.0);
        }
        std::cout << "FFT...\n";
        waves = FFT(waves);
        std::cout << "truncating smallest...\n";
        std::vector<std::tuple<std::complex<double>, size_t>> filterlist(waves.size());

        for (size_t i = 0; i < waves.size(); i++) {
            filterlist[i] = std::make_tuple(waves[i], i);
        }

        std::sort(filterlist.begin(), filterlist.end(), [](std::tuple<std::complex<double>, size_t> a, std::tuple<std::complex<double>, size_t> b) {return abs(std::get<0>(a)) > abs(std::get<0>(b));});
        
        for (size_t i = (double)waves.size() / keepfrac; i < waves.size(); i++) {
            std::get<0>(filterlist[i]) = 0;
        }

        std::sort(filterlist.begin(), filterlist.end(), [](std::tuple<std::complex<double>, size_t> a, std::tuple<std::complex<double>, size_t> b) {return std::get<1>(a) < std::get<1>(b);});

        for (size_t i = 0; i < waves.size(); i++) {
            waves[i] = std::get<0>(filterlist[i]);
        }

        std::cout << "reconstructing...\n";
        waves = IFFT(waves);

        std::cout << "decoding...\n";

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int index = 3 * (j + width * i);
                int hilbidx = 3 * gilbidx(j, i, width, height);

                for (int k = 0; k < channels; k++) {
                    rawOut[index + k] = std::min(255, (int)abs(waves[hilbidx + k]));
                }
            }
        }

        std::cout << "writing...\n";
        std::string outpath = "out" + std::to_string((long)keepfrac) + ".ppm";
        savePPM(outpath, width, height, rawOut);
        keepfrac *= 2;
    }
}
