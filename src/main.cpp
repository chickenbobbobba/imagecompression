#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>
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
        unsigned char* temp = stbi_load(path.c_str(), &width, &height, &channels, 3);
        if (temp) {
            image2d.assign(temp, temp + width * height * channels);
        } else {
            std::cout << "error loading image!\n";
            throw;
        }
        std::cout << "width: " << width << "\nheight: " << height << "\nchannels: " << channels << "\n";
    }

    void savePPM(const std::string& path) {
        assert(image2d.size() == width * height * channels);
        std::ofstream out(path, std::ios::binary);
        out << "P6\n" << width << " " << height << "\n255\n";
        out.write(reinterpret_cast<char*>(image2d.data()), image2d.size());
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

        for (char i : image1d) {
            waves.emplace_back((double)i, 0.0);
        }

        waves = FFT(waves);
        for (auto& i : waves) {
            std::cout << i << "\n";
        }
    }

    void wavesToImage1d() {
        if (image1d.size() < waves.size())
            image1d.resize(waves.size());

        waves = IFFT(waves);

        for (size_t i = 0; i < waves.size(); i++) {
            image1d[i] = abs(waves[i]);
        }
    }
};

int main(int argc, char** argv){
    double keepfrac = 1.0f;

    Image image;
    std::string filepath = "../resources/mandelbrot.png";
    if (argc > 1) filepath = argv[1];
    std::cout << "filepath: " << filepath << "\n";
    std::cout << "reading...\n";
    image.loadImage(filepath);
    std::cout << "encoding...\n";
    image.hilbTo1d();
    std::cout << "FFT...\n";
    image.image1dToWaves();
    std::cout << "truncating smallest...\n";
    // {

    //     std::vector<std::tuple<std::complex<double>, size_t>> filterlist(image.image1d.size());
    //     for (size_t i = 0; i < filterlist.size(); i++) {
    //         filterlist[i] = std::make_tuple(image.image1d[i], i);
    //     }
    //     std::sort(filterlist.begin(), filterlist.end(), [](const std::tuple<std::complex<double>&, size_t>& a, const std::tuple<std::complex<double>&, size_t>& b) {return abs(std::get<0>(a)) > abs(std::get<0>(b));});
    //     for (size_t i = (double)filterlist.size() / keepfrac; i < filterlist.size(); i++) {
    //         std::get<0>(filterlist[i]) = 0;
    //     }
        
    //     std::sort(filterlist.begin(), filterlist.end(), [](const std::tuple<std::complex<double>&, size_t>& a, const std::tuple<std::complex<double>&, size_t>& b) {return std::get<1>(a) < std::get<1>(b);});
    //     for (size_t i = 0; i < filterlist.size(); i++) {
    //         image.image1d[i] = std::min(255, (int)std::get<0>(filterlist[i]).real());
    //     }
    
    // }
    std::cout << "reconstructing...\n";
    image.wavesToImage1d();

    std::cout << "decoding...\n";
    image.hilbTo2d();
    std::cout << "writing...\n";
    std::string outpath = "out.ppm";
    image.savePPM(outpath);

}
