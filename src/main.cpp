#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>
#include <hilbert.hpp>
#include <fstream>

#include <fftw3.h>
#include <fftwrap.hpp>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include <threadpool.h>

class Subsect {
    private:
    public:

    std::vector<unsigned char> Y;
    std::vector<unsigned char> Cb;
    std::vector<unsigned char> Cr;

    std::vector<unsigned char> raw;
    std::vector<std::complex<double>> waves;

    size_t start;
    size_t end;

    void assignRawData(const std::vector<unsigned char>& data, size_t start, size_t endExclusive) {
        if (start > endExclusive || endExclusive > data.size()) {
            throw std::out_of_range("assignRawData: invalid range");
        }
        raw.assign(data.begin() + start, data.begin() + endExclusive);
        FFT::init(raw.size());
        this->start = start;
        this->end = endExclusive;
    }
    
    void integrateRawData(std::vector<unsigned char>& data) {
        assert(data.size() >= raw.size() + start);
        memcpy(data.data() + start, raw.data(), raw.size());
    }
    
    std::vector<std::complex<double>> toWaves(const std::vector<unsigned char>& data) {
        std::vector<std::complex<double>> waves;
        waves.reserve(data.size());
        for (unsigned char v : data) {
            waves.emplace_back((double)v, 0.0);
        }

        FFT::forward(waves);
        return waves; 
    }

    std::vector<unsigned char> toData(const std::vector<std::complex<double>>& data) {
        std::vector<unsigned char> raw(data.size());

        auto temp = data;
        FFT::backward(temp);

        for (size_t i = 0; i < temp.size(); i++) {
            raw[i] = (unsigned char) std::clamp((long)(abs(temp[i])/data.size()), 0L, 255L);
        }

        return raw;
    }
};

class Image {
public:
    int width;
    int height;
    int channels;
    long segLen;
    size_t length;
    size_t rawLength;

    std::vector<unsigned char> rawData; // standard linear mapping
    std::vector<unsigned char> hilbMap; // hilbert mapping

    std::vector<Subsect> subsects;

    void loadImage(const std::string& path) {
        unsigned char* temp = stbi_load(path.c_str(), &width, &height, &channels, 3);
        if (temp) {
            channels = 3;
            length = width * height;
            rawLength = length * channels;
            rawData.resize(rawLength);
            rawData.assign(temp, temp + rawLength);
            STBI_FREE(temp);
        } else {
            std::cout << "error loading image!\n";
            throw;
        }
        std::cout << "width: " << width << "\nheight: " << height << "\npixels: " << length << "\nchannels: " << channels << "\n";
    }

    void subdivide(std::vector<unsigned char> data, long count) {
        for (long i = 0; i < count; i++) {
            Subsect temp;
            long start = i * data.size()/count;
            long endExclusive = (i+1) * data.size()/count; // exclusive end
            temp.assignRawData(data, start, endExclusive);
            // record the segment bounds
            temp.start = start;
            temp.end = endExclusive;
            subsects.push_back(std::move(temp));
        }
    }

    void savePPM(const std::string& path) {
        std::ofstream out(path, std::ios::binary);
        out << "P6\n" << width << " " << height << "\n255\n";
        out.write(reinterpret_cast<char*>(rawData.data()), rawData.size());
    }

    void rawToHilb() {
        hilbMap.resize(rawLength);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int index = channels * (j + width * i);
                int hilbidx = channels * gilbidx(j, i, width, height);

                for (int k = 0; k < channels; k++) {
                    hilbMap[hilbidx + k] = rawData[index + k];
                }
            }
        }
    }

    void hilbToRaw() {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int index = channels * (j + width * i);
                int hilbidx = channels * gilbidx(j, i, width, height);

                for (int k = 0; k < channels; k++) {
                    rawData[index + k] = hilbMap[hilbidx + k];
                }
            }
        }
    }
    
};

int main(int argc, char** argv) {
    ThreadPool pool(std::thread::hardware_concurrency());
    
    Image image;
    std::string filepath = "";
    if (argc > 1) filepath = argv[1];
    std::cout << "filepath: " << filepath << "\n";
    image.loadImage(filepath);
    image.rawToHilb();
    image.subdivide(image.hilbMap, (long)sqrt(image.rawLength));

    for (auto &i : image.subsects) {
        i.waves = i.toWaves(i.raw);

        double avg = 0; // RMS average
        for (auto i : i.waves) avg += pow(abs(i), 2);
        avg /= i.waves.size();
        avg = sqrt(avg);

        for (auto& i : i.waves) i *= (abs(i) > avg);

        i.raw = i.toData(i.waves);
        i.integrateRawData(image.hilbMap);
    }
    image.hilbToRaw();
    image.savePPM("img.ppm");
}
