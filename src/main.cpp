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

class Image {
public:
    int width;
    int height;
    int channels;
    std::vector<unsigned char> image2d; // standard linear mapping
    std::vector<unsigned char> image1d; // hilbert mapping
    std::vector<std::complex<double>> waves;

    void loadImage(const std::string& path) {
        image2d = stbi_load(path.c_str(), &width, &height, &channels, 3);
        std::cout << "width: " << width << "\nheight: " << height << "\nchannels: " << channels << "\n";
    }

    void savePPM(const std::string& path) {
        assert(image2d.size() == width * height * channels);
        std::ofstream out(path, std::ios::binary);
        out << "P6\n" << width << " " << height << " " << channels << "\n255\n";
        out.write((const char*)image2d.data(), image2d.size());
    }

    void hilbTo1d() {
        image1d.resize(width * height * channels);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int index = channels * (j + width * i);
                int hilbidx = channels * gilbidx(j, i, width, height);

                for (int k = 0; k < channels; k++) {
                    image1d[hilbidx + k] = image2d[index + k];
                }
            }
        }
    }

    void hilbTo2d() {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int index = channels * (j + width * i);
                int hilbidx = channels * gilbidx(j, i, width, height);

                for (int k = 0; k < channels; k++) {
                    image2d[index + k] = image1d[hilbidx + k];
                }
            }
        }
    }

    void image1dToWaves() {
        waves.reserve(image1d.size());

        for (auto i : image1d) {
            waves.emplace_back(i, 0.0);
        }

        waves = FFT(image1d);
    }

    void wavesToImage1d() {
        if (image1d.size() < waves.size())
            image1d.resize(waves.size());

        image1d = IFFT(waves);
    }
};

int main(int argc, char** argv){
    double keepfrac = /* 1/ */1;
    while (keepfrac <= 1) {
        Image image;
        std::string filepath = "../resources/mandelbrot.png";
        if (argc > 1) filepath = argv[1];
        std::cout << "filepath: " << filepath << "\n";
        std::cout << "reading...\n";
        std::cout << "encoding...\n";
        std::cout << "FFT prep...\n";
        std::cout << "FFT...\n";
        std::cout << "truncating smallest...\n";
        std::vector<std::tuple<std::complex<double>, size_t>> filterlist(waves.size());

        for (size_t i = 0; i < waves.size(); i++) {
            filterlist[i] = std::make_tuple(waves[i], i);
        }

        std::sort(filterlist.begin(), filterlist.end(), [](const std::tuple<std::complex<double>&, size_t>& a, const std::tuple<std::complex<double>&, size_t>& b) {return abs(std::get<0>(a)) > abs(std::get<0>(b));});
        
        for (size_t i = (double)waves.size() / keepfrac; i < waves.size(); i++) {
            std::get<0>(filterlist[i]) = 0;
        }

        std::cout << "reconstructing...\n";

        std::sort(filterlist.begin(), filterlist.end(), [](const std::tuple<std::complex<double>&, size_t>& a, const std::tuple<std::complex<double>&, size_t>& b) {return std::get<1>(a) < std::get<1>(b);});

        for (size_t i = 0; i < waves.size(); i++) {
            waves[i] = std::get<0>(filterlist[i]);
        }

        std::cout << "decoding...\n";
        std::cout << "writing...\n";
        std::string outpath = "out" + std::to_string((long)keepfrac) + ".ppm";
        keepfrac *= 4;
    }
}
