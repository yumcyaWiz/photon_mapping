#include "camera.h"
#include "image.h"
#include "primitive.h"

int main() {
  const int width = 512;
  const int height = 512;

  const auto shape = std::make_shared<Sphere>(Vec3(0), 1.0f);
  const Primitive sphere(shape, nullptr);

  Camera camera(Vec3(0, 0, 3), Vec3(0, 0, -1));

  Image image(width, height);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      const float u = (2.0f * j - width) / height;
      const float v = (2.0f * i - height) / height;

      Ray ray;
      float pdf;
      if (camera.sampleRay(u, v, ray, pdf)) {
        IntersectInfo info;
        if (sphere.intersect(ray, info)) {
          image.setPixel(i, j, 0.5f * (info.hitNormal + 1.0f));
        }
      } else {
        image.setPixel(i, j, Vec3(0));
      }
    }
  }

  image.writePPM("output.ppm");

  return 0;
}