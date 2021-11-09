#include "camera.h"
#include "image.h"
#include "primitive.h"
#include "scene.h"

int main() {
  const int width = 512;
  const int height = 512;

  const auto sphere_shape = std::make_shared<Sphere>(Vec3(0, 1, 0), 1.0f);
  const auto sphere2_shape = std::make_shared<Sphere>(Vec3(-1, 1, -1), 1.0f);
  const auto sphere3_shape = std::make_shared<Sphere>(Vec3(1, 1, 1), 1.0f);
  const auto floor_shape =
      std::make_shared<Plane>(Vec3(5, 0, -5), Vec3(-10, 0, 0), Vec3(0, 0, 10));

  const Primitive floor(floor_shape, nullptr);
  const Primitive sphere(sphere_shape, nullptr);
  const Primitive sphere2(sphere2_shape, nullptr);
  const Primitive sphere3(sphere3_shape, nullptr);

  Scene scene;
  scene.addPrimitive(floor);
  scene.addPrimitive(sphere);
  scene.addPrimitive(sphere2);
  scene.addPrimitive(sphere3);
  scene.build();

  Camera camera(Vec3(0, 1, 5), Vec3(0, 0, -1));

  Image image(width, height);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      const float u = (2.0f * j - width) / height;
      const float v = (2.0f * i - height) / height;

      Ray ray;
      float pdf;
      if (camera.sampleRay(Vec2(u, v), ray, pdf)) {
        IntersectInfo info;
        if (scene.intersect(ray, info)) {
          image.setPixel(i, j, 0.5f * (info.surfaceInfo.normal + 1.0f));
        }
      } else {
        image.setPixel(i, j, Vec3(0));
      }
    }
  }

  image.writePPM("output.ppm");

  return 0;
}