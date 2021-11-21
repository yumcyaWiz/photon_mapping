#include "camera.h"
#include "image.h"
#include "scene.h"

int main() {
  const int width = 512;
  const int height = 512;

  Scene scene;
  scene.loadModel("CornellBox-Original.obj");
  scene.build();

  Camera camera(Vec3f(0, 1, 7), Vec3f(0, 0, -1), 0.25f * PI);

  Image image(width, height);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      const float u = (2.0f * j - width) / height;
      const float v = (2.0f * i - height) / height;

      Ray ray;
      float pdf;
      if (camera.sampleRay(Vec2f(u, v), ray, pdf)) {
        IntersectInfo info;
        if (scene.intersect(ray, info)) {
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