#ifndef _IMAGE_H
#define _IMAGE_H
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "core.h"

class Image {
 private:
  unsigned int width;
  unsigned int height;
  std::vector<float> pixels;

 public:
  Image(unsigned int width, unsigned int height)
      : width(width), height(height) {
    pixels.resize(3 * width * height);
  }

  Vec3 getPixel(unsigned int i, unsigned int j) const {
    const unsigned int idx = 3 * j + 3 * width * i;
    return Vec3(pixels[idx], pixels[idx + 1], pixels[idx + 2]);
  }

  void setPixel(unsigned int i, unsigned int j, const Vec3& rgb) {
    const unsigned int idx = 3 * j + 3 * width * i;
    pixels[idx] = rgb[0];
    pixels[idx + 1] = rgb[1];
    pixels[idx + 2] = rgb[2];
  }

  void writePPM(const std::string& filename) {
    std::ofstream file(filename);

    file << "P3" << std::endl;
    file << width << " " << height << std::endl;
    file << "255" << std::endl;

    for (unsigned int i = 0; i < height; ++i) {
      for (unsigned int j = 0; j < width; ++j) {
        const Vec3 rgb = getPixel(i, j);
        const unsigned int R =
            std::clamp(static_cast<unsigned int>(255.0f * rgb[0]), 0u, 255u);
        const unsigned int G =
            std::clamp(static_cast<unsigned int>(255.0f * rgb[1]), 0u, 255u);
        const unsigned int B =
            std::clamp(static_cast<unsigned int>(255.0f * rgb[2]), 0u, 255u);
        file << R << " " << G << " " << B << std::endl;
      }
    }

    file.close();
  }
};

#endif