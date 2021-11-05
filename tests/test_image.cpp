#include "image.h"

int main() {
  const int width = 512;
  const int height = 512;

  Image image(width, height);
  for (int i = 0; i < height; ++i) {
    const float u = static_cast<float>(i) / height;
    for (int j = 0; j < width; ++j) {
      const float v = static_cast<float>(j) / height;
      image.setPixel(i, j, Vec3(u, v, 1.0f));
    }
  }

  image.writePPM("output.ppm");

  return 0;
}