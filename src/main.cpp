#include <algorithm>
#include <complex>
#include <iostream>
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

    // std::sort(waves.begin(), waves.end(), [](std::complex<double> a, std::complex<double> b){return abs(a) > abs(b);});
    // for (size_t i = waves.size() / 10; i < waves.size(); i++) {
    //     waves[i] = 0;
    // }

    std::cout << "reconstructing...\n";
    waves = IFFT(waves);

    std::cout << "decoding...\n";

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int index = 3 * (j + width * i);
            int hilbidx = 3 * gilbidx(j, i, width, height);

            for (int k = 0; k < channels; k++) {
                rawOut[index + k] = abs(waves[hilbidx + k]);
            }
        }
    }

    std::cout << "writing...\n";
    std::string outpath = "out.ppm";
    savePPM(outpath, width, height, rawOut);
}
