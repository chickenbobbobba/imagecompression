#include <algorithm>
#include <climits>
#include <cmath>
#include <complex>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <hilbert.hpp>
#include <fstream>

#include <fftw3.h>
#include <fftwrap.hpp>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>



class Image {
public:
    int width;
    int height;
    int channels = 3;
    long segLen;
    long colreduc;
    size_t length;
    size_t rawLength;

    std::vector<unsigned char> rawData; // standard linear mapping
    std::vector<unsigned char> image1d; // hilbert mapping

    class Segment {
        private:
        Image& parent;
        public:
        Segment(Image& parent_, size_t segLen, size_t index, long yScale, long cScale) : parent(parent_){
            Y.resize(segLen/yScale);
            Cb.resize(segLen/cScale);
            Cr.resize(segLen/cScale);

            if (parent.image1d.size() == 0) throw;

            //fill Y
            for (auto i = index; i < index + segLen; i += parent.channels * yScale) {
                Y.push_back(0.299 * parent.image1d[i] + 0.587 * parent.image1d[i+1] + 0.114 * parent.image1d[i+2]);
            }
            //fill Cb Cr
            for (auto i = index; i < index + segLen; i += parent.channels * cScale) {
                Cb.push_back(128 - 0.168736 * parent.image1d[i] - 0.331264 * parent.image1d[i+1] + 0.5 * parent.image1d[i+2]);
                Cr.push_back(128 + 0.5 * parent.image1d[i] - 0.418688 * parent.image1d[i+1] - 0.081312 * parent.image1d[i+2]);
            }
        }
        std::vector<unsigned char> Y;
        std::vector<unsigned char> Cb;
        std::vector<unsigned char> Cr;

        std::vector<std::complex<double>> wavY;
        std::vector<std::complex<double>> wavCb;
        std::vector<std::complex<double>> wavCr;

        // std::vector<std::complex<double>> getFreqs(const std::vector<unsigned char>& data) {
        //     std::vector<std::complex<double>> waves(data.size());
        //     for (unsigned char i : data) {
        //         waves.emplace_back((double)i, 0.0);
        //     }
        //     return FFT::forward(waves);
        // }

        // std::vector<unsigned char> getData(const std::vector<std::complex<double>>& data) {
        //     std::vector<unsigned char> raw(data.size());

        //     auto temp = FFT::backward(data);

        //     for (size_t i = 0; i < temp.size(); i++) {
        //         raw[i] = (unsigned char)std::min(255L, (long)abs(temp[i]));
        //     }
        //     return raw;
        // }

        // std::vector<std::complex<double>> compress(const std::vector<unsigned char>& data) {
        //     double twiddle = 1;

        //     auto waves = getFreqs(data);

        //     double avg = 0;
        //     for (const auto& i : waves) {
        //         avg += abs(i);
        //     }
        //     avg /= waves.size() * twiddle;

        //     for (auto& i : waves) {
        //         if (abs(i) < avg) i *= 0;
        //     }
        //     return waves;
        // }

        // void compressAll() {
        //     wavY = compress(Y);
        //     wavCb = compress(Cb);
        //     wavCr = compress(Cr);
        // }
    };

    std::vector<Segment> segments;

    // void compressSegments() {
    //     for (auto& i : segments) {
    //         i.compressAll();
    //     }
    // }
    
    void loadImage(const std::string& path) {
        unsigned char* temp = stbi_load(path.c_str(), &width, &height, &channels, 3);
        if (temp) {
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

    void segmentRaw() {
        segLen = 64 * channels * sqrt(length); //todo : make adaptive probably
        size_t i = 0;
        long count = 0;
        while (i < rawLength) {
            count++;
            if (i + segLen > rawLength) segLen = rawLength - i;
            segments.emplace_back(*this, segLen, i, 1, 16);
            i += segLen;
        }
        std::cout << "created " << count << " segments\n";
    }

    void rejoinSegments() {

    }

    void saveNift(const std::string& path) {
        /*
        spec:
        Nxn sized hIlbert and Fourier Transform
        greyscale, rgb and rgba support
        per image variable size transformed coeffs

        bits : what its for

        u4 : block size is (image size)^(1/2)
        u32 : width
        u32 : height
        u2 : format (0 = grey, 2 = rgb, 3 = rgba)

        u4 : num bits for transform real coeffs -1 (so 1111 -> 10000 -> 16 bits of precision, 0000 -> 0001 -> 1 bit of precision)
        u4 : same for imag coeffs

        */
    }

    void savePPM(const std::string& path) {
        std::ofstream out(path, std::ios::binary);
        out << "P6\n" << width << " " << height << "\n255\n";
        out.write(reinterpret_cast<char*>(rawData.data()), rawData.size());
    }

    void hilbTo1d() {
        image1d.resize(rawLength);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int index = channels * (j + width * i);
                int hilbidx = channels * gilbidx(j, i, width, height);

                for (int k = 0; k < channels; k++) {
                    image1d[hilbidx + k] = rawData[index + k];
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
                    rawData[index + k] = image1d[hilbidx + k];
                }
            }
        }
    }

};

int main(int argc, char** argv) {
    
    Image image;
    std::string filepath = "";
    if (argc > 1) filepath = argv[1];
    std::cout << "filepath: " << filepath << "\n";
    image.loadImage(filepath);

    FFT::init((size_t)image.rawData.size());
    
    std::vector<std::complex<double>> b;
    std::vector<unsigned char> a(image.rawData.size());
    for (int i = 0; i < image.height; i++) {
        for (int j = 0; j < image.width; j++) {
            int index = image.channels * (j + image.width * i);
            int hilbidx = image.channels * gilbidx(j, i, image.width, image.height);

            for (int k = 0; k < image.channels; k++) {
                a[hilbidx + k] = image.rawData[index + k];
            }
        }
    }
    b.reserve(image.rawData.size());
    for (unsigned char i : a) {
        b.emplace_back((double)i, 0.0);
    }

    
    FFT::forward(b);
    auto c = b;

    long frames = 20*60;
    for (int i = 0; i < frames; i++) {
        b = c;
        std::cout << b.size() - b.size() * i / frames << "\n";
        for (int j = b.size() - b.size() * i / frames; j < b.size(); j++) {
            b[j] = 0;
        }
        c = b;

        FFT::backward(b);

        for (int i = 0; i < image.height; i++) {
            for (int j = 0; j < image.width; j++) {
                int index = image.channels * (j + image.width * i);
                int hilbidx = image.channels * gilbidx(j, i, image.width, image.height);

                for (int k = 0; k < image.channels; k++) {
                    double val = std::abs(b[hilbidx + k])/(double)b.size();
                    image.rawData[index + k] = (unsigned char)std::clamp(val, 0.0, 255.0);
                }
            }
        }

        std::string num = std::to_string(i);
        while (num.size() < 5) num = "0" + num;
        image.savePPM("frame_" + num + ".ppm");
    }
}
