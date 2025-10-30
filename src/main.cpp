#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <thread>
#include <unistd.h>
#include <utility>
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
    std::vector<std::complex<double>> CbCr; // joined to save FFT operations

    std::vector<unsigned char> raw;

    size_t start;
    size_t end;

    void assignRawData(const std::vector<unsigned char>& data, size_t start, size_t end) {
        if (start > end || end > data.size()) {
            throw std::out_of_range("assignRawData: invalid range");
        }
        raw.assign(data.begin() + start, data.begin() + end);
        FFT::init(raw.size());
        this->start = start;
        this->end = end;
    }
    
    void integrateRawData(std::vector<unsigned char>& data) {
        assert(data.size() >= raw.size() + start);
        memcpy(data.data() + start, raw.data(), raw.size());
    }

    void toYCbCr(size_t yScale, size_t cScale) {
        size_t segLen = raw.size();
        size_t numPixels = segLen / 3;
        
        Y.reserve(numPixels / yScale);
        CbCr.reserve(numPixels / cScale);
        
        // fill y
        for (size_t px = 0; px < numPixels; px += yScale) {
            double r_sum = 0, g_sum = 0, b_sum = 0;
            size_t count = 0;
            
            for (size_t j = 0; j < yScale && (px + j) < numPixels; j++) {
                size_t idx = (px + j) * 3;
                r_sum += pow(raw[idx + 0], 2.2);
                g_sum += pow(raw[idx + 1], 2.2);
                b_sum += pow(raw[idx + 2], 2.2);
                count++;
            }
            
            double r_avg = pow(r_sum / count, 1.0/2.2);
            double g_avg = pow(g_sum / count, 1.0/2.2);
            double b_avg = pow(b_sum / count, 1.0/2.2);
            
            Y.push_back(0.299 * r_avg + 0.587 * g_avg + 0.114 * b_avg);
        }
        
        // fill Cb Cr
        for (size_t px = 0; px < numPixels; px += cScale) {
            double r_sum = 0, g_sum = 0, b_sum = 0;
            size_t count = 0;
            
            for (size_t j = 0; j < cScale && (px + j) < numPixels; j++) {
                size_t idx = (px + j) * 3;
                r_sum += pow(raw[idx + 0], 2.2);
                g_sum += pow(raw[idx + 1], 2.2);
                b_sum += pow(raw[idx + 2], 2.2);
                count++;
            }
            
            double r_avg = pow(r_sum / count, 1.0/2.2);
            double g_avg = pow(g_sum / count, 1.0/2.2);
            double b_avg = pow(b_sum / count, 1.0/2.2);
            
            CbCr.emplace_back(128 - 0.168736 * r_avg - 0.331264 * g_avg + 0.5 * b_avg,
                              128 + 0.5 * r_avg - 0.418688 * g_avg - 0.081312 * b_avg);
        }
    }

    void fromYCbCr() {
        size_t segLen = end - start;
        if (segLen == 0) return;
        size_t np = segLen / 3;

        if (Y.empty() && CbCr.empty()) {
            raw.assign(segLen, 0);
            return;
        }

        // upscale CbCr if needed
        if (!CbCr.empty() && CbCr.size() != np) {
            size_t oldSize = CbCr.size();
            size_t half = (oldSize + 1) / 2;
            
            FFT::init(oldSize);
            FFT::forward(CbCr);
            
            // normalize after forward FFT
            for (auto& f : CbCr) f /= (double)oldSize;
            
            // create new padded array
            std::vector<std::complex<double>> temp(np, {0.0, 0.0});
            
            // copy low (positive) frequencies to beginning
            for (size_t i = 0; i < half; i++) {
                temp[i] = CbCr[i];
            }
            
            // copy high (negative) frequencies to end
            for (size_t i = half; i < oldSize; i++) {
                temp[np - (oldSize - i)] = CbCr[i];
            }
            
            FFT::init(np);
            FFT::backward(temp);
            
            CbCr = std::move(temp);
        }

        // upscale Y if needed
        if (!Y.empty() && Y.size() != np) {
            size_t oldY = Y.size();
            size_t halfY = (oldY + 1) / 2;
            auto tempY = toWaves(Y);
            size_t NY = tempY.size();
            
            std::vector<std::complex<double>> paddedY(np, {0.0, 0.0});
            
            for (size_t i = 0; i < halfY; i++) {
                paddedY[i] = tempY[i];
            }

            for (size_t i = halfY; i < NY; i++) {
                paddedY[np - (NY - i)] = tempY[i];
            }

            FFT::init(np);
            Y = toData(paddedY);
        }

        raw.assign(segLen, 0);

        for (size_t px = 0; px < np; ++px) {
            size_t y_idx = (Y.empty() ? 0 : std::min(px, Y.size() - 1));
            size_t c_idx = (CbCr.empty() ? 0 : std::min(px, CbCr.size() - 1));
            
            double Yv = (Y.empty() ? 0.0 : (double)Y[y_idx]);
            double Cbv = (CbCr.empty() ? 128.0 : CbCr[c_idx].real());
            double Crv = (CbCr.empty() ? 128.0 : CbCr[c_idx].imag());
            
            double R = Yv + 1.402   * (Crv - 128.0);
            double G = Yv - 0.344136* (Cbv - 128.0) - 0.714136 * (Crv - 128.0);
            double B = Yv + 1.772   * (Cbv - 128.0);

            // std::cout << R << " " << G << " " << B << "\n";
            
            auto clamp_byte = [](double v) -> unsigned char {
                return static_cast<unsigned char>(std::clamp((int)std::lround(v), 0, 255));
            };
            
            size_t base = px * 3;
            raw[base + 0] = clamp_byte(R);
            raw[base + 1] = clamp_byte(G);
            raw[base + 2] = clamp_byte(B);
        }
    }
    
    std::vector<std::complex<double>> toWaves(const std::vector<unsigned char>& data) {
        std::vector<std::complex<double>> waves;
        waves.reserve(data.size());
        for (unsigned char v : data) {
            waves.emplace_back((double)v, 0.0);
        }
        
        FFT::init(waves.size());
        FFT::forward(waves);
        for (auto& i : waves) i /= (double)data.size();
        return waves; 
    }

    std::vector<unsigned char> toData(const std::vector<std::complex<double>>& data) {
        std::vector<unsigned char> raw(data.size());

        auto temp = data;
        FFT::init(temp.size());
        FFT::backward(temp);

        for (size_t i = 0; i < temp.size(); i++) {
            raw[i] = (unsigned char) std::clamp((long)(abs(temp[i])), 0L, 255L);
        }

        return raw;
    }

    std::vector<unsigned char> encodeBlockData() {
        /*
        32 : width
        32 : height
        32 : index of 2x8 bit
        32 : index of 2x6 bit
        32 : index of 2x4 bit
        32 : index of 2x2 bit
        */
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
            start -= start % 3;
            long endExclusive = (i+1) * data.size()/count;
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

    char mode = argv[2][0];

    Image image;
    std::string filepath = "";
    if (argc > 1) filepath = argv[1];
    std::cout << "filepath: " << filepath << "\n";
    image.loadImage(filepath);
    
    // encode
    image.rawToHilb();
    image.subdivide(image.hilbMap, sqrt(image.rawLength) * 16);

    for (auto &i : image.subsects) {
        i.toYCbCr(1, 4);
        
        // auto wavesY = i.toWaves(i.Y);
        // FFT::init(i.CbCr.size());
        // FFT::forward(i.CbCr);
        // for (auto& a : i.CbCr) a /= i.CbCr.size();

        // i.Y = i.toData(wavesY);
        // FFT::backward(i.CbCr);
        
        if (mode == 'c') {
            for (auto& a : i.Y) a = 128.0;
        }
        if (mode == 'y') {
            for (auto& a : i.CbCr) a = std::complex<double>(128, 128);
        }

        i.fromYCbCr();
        i.integrateRawData(image.hilbMap);
    }
    
    image.hilbToRaw();
    image.savePPM("img.ppm");
}
