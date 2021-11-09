#include "camera.h"
#include "image.h"

int main() {
  const int width = 512;
  const int height = 512;

  Camera camera(Vec3(0), Vec3(0, 0, -1));

  Image image(width, height);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      const float u = (2.0f * j - width) / height;
      const float v = (2.0f * i - height) / height;

      Ray ray;
      float pdf;
      if (camera.sampleRay(Vec2(u, v), ray, pdf)) {
        image.setPixel(i, j, 0.5f * (ray.direction + 1.0f));
      } else {
        image.setPixel(i, j, Vec3(0));
      }
    }
  }

  image.writePPM("output.ppm");

  return 0;
}