#include "camera.h"
#include "image.h"
#include "primitive.h"

int main() {
  const int width = 512;
  const int height = 512;
  Image image(width, height);

  const auto shape = std::make_shared<Plane>(Vec3f(2.5, 0, -2.5),
                                             Vec3f(-5, 0, 0), Vec3f(0, 0, 5));
  const Primitive plane(shape, nullptr);

  Camera camera(Vec3f(0, 3, 3), Vec3f(0, 0, -1));

  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      const float u = (2.0f * j - width) / height;
      const float v = (2.0f * i - height) / height;

      Ray ray;
      float pdf;
      if (camera.sampleRay(Vec2f(u, v), ray, pdf)) {
        IntersectInfo info;
        if (plane.intersect(ray, info)) {
          image.setPixel(i, j, 0.5f * (info.surfaceInfo.normal + 1.0f));
        }
      } else {
        image.setPixel(i, j, Vec3f(0));
      }
    }
  }

  image.writePPM("output.ppm");

  return 0;
}