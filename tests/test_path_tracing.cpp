#include <omp.h>

#include "camera.h"
#include "image.h"
#include "integrator.h"
#include "photon_map.h"
#include "primitive.h"
#include "scene.h"

int main() {
  const int width = 512;
  const int height = 512;
  const int n_samples = 10000;
  const int max_depth = 100;
  const Vec3 camPos(2.78, 2.73, -9);
  const Vec3 lookAt(2.78, 2.73, 2.796);

  Image image(width, height);

  const Camera camera(camPos, normalize(lookAt - camPos), 0.25 * PI);

  // cornellbox scene
  const auto white = std::make_shared<Lambert>(Vec3(0.8));
  const auto red = std::make_shared<Lambert>(Vec3(0.8, 0.05, 0.05));
  const auto green = std::make_shared<Lambert>(Vec3(0.05, 0.8, 0.05));
  const auto mirror = std::make_shared<Mirror>(Vec3(0.99));

  const auto floor =
      std::make_shared<Plane>(Vec3(0), Vec3(0, 0, 5.592), Vec3(5.56, 0, 0));
  const auto rightWall =
      std::make_shared<Plane>(Vec3(0), Vec3(0, 5.488, 0), Vec3(0, 0, 5.592));
  const auto leftWall = std::make_shared<Plane>(
      Vec3(5.56, 0, 0), Vec3(0, 0, 5.592), Vec3(0, 5.488, 0));
  const auto ceil = std::make_shared<Plane>(Vec3(0, 5.488, 0), Vec3(5.56, 0, 0),
                                            Vec3(0, 0, 5.592));
  const auto backWall = std::make_shared<Plane>(
      Vec3(0, 0, 5.592), Vec3(0, 5.488, 0), Vec3(5.56, 0, 0));

  const auto shortBox1 = std::make_shared<Plane>(
      Vec3(1.3, 1.65, 0.65), Vec3(-0.48, 0, 1.6), Vec3(1.6, 0, 0.49));
  const auto shortBox2 = std::make_shared<Plane>(
      Vec3(2.9, 0, 1.14), Vec3(0, 1.65, 0), Vec3(-0.5, 0, 1.58));
  const auto shortBox3 = std::make_shared<Plane>(
      Vec3(1.3, 0, 0.65), Vec3(0, 1.65, 0), Vec3(1.6, 0, 0.49));
  const auto shortBox4 = std::make_shared<Plane>(
      Vec3(0.82, 0, 2.25), Vec3(0, 1.65, 0), Vec3(0.48, 0, -1.6));
  const auto shortBox5 = std::make_shared<Plane>(
      Vec3(2.4, 0, 2.72), Vec3(0, 1.65, 0), Vec3(-1.58, 0, -0.47));

  const auto tallBox1 = std::make_shared<Plane>(
      Vec3(4.23, 3.30, 2.47), Vec3(-1.58, 0, 0.49), Vec3(0.49, 0, 1.59));
  const auto tallBox2 = std::make_shared<Plane>(
      Vec3(4.23, 0, 2.47), Vec3(0, 3.3, 0), Vec3(0.49, 0, 1.59));
  const auto tallBox3 = std::make_shared<Plane>(
      Vec3(4.72, 0, 4.06), Vec3(0, 3.3, 0), Vec3(-1.58, 0, 0.5));
  const auto tallBox4 = std::make_shared<Plane>(
      Vec3(3.14, 0, 4.56), Vec3(0, 3.3, 0), Vec3(-0.49, 0, -1.6));
  const auto tallBox5 = std::make_shared<Plane>(
      Vec3(2.65, 0, 2.96), Vec3(0, 3.3, 0), Vec3(1.58, 0, -0.49));

  const auto light_shape = std::make_shared<Plane>(
      Vec3(3.43, 5.486, 2.27), Vec3(0, 0, 1.05), Vec3(-1.3, 0, 0));
  const auto light = std::make_shared<AreaLight>(Vec3(34, 19, 10), light_shape);

  Scene scene;
  scene.addPrimitive(Primitive(floor, mirror));
  scene.addPrimitive(Primitive(rightWall, red));
  scene.addPrimitive(Primitive(leftWall, green));
  scene.addPrimitive(Primitive(ceil, white));
  scene.addPrimitive(Primitive(backWall, white));
  scene.addPrimitive(Primitive(shortBox1, white));
  scene.addPrimitive(Primitive(shortBox2, white));
  scene.addPrimitive(Primitive(shortBox3, white));
  scene.addPrimitive(Primitive(shortBox4, white));
  scene.addPrimitive(Primitive(shortBox5, white));
  scene.addPrimitive(Primitive(tallBox1, white));
  scene.addPrimitive(Primitive(tallBox2, white));
  scene.addPrimitive(Primitive(tallBox3, white));
  scene.addPrimitive(Primitive(tallBox4, white));
  scene.addPrimitive(Primitive(tallBox5, white));
  scene.addPrimitive(Primitive(light_shape, white, light));
  scene.build();

  // photon tracing and build photon map
  PathTracing integrator(max_depth);

#pragma omp parallel for collapse(2) schedule(dynamic, 1)
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      // init sampler
      UniformSampler sampler(j + width * i);

      for (int k = 0; k < n_samples; ++k) {
        const float u = (2.0f * (j + sampler.getNext1D()) - width) / height;
        const float v = (2.0f * (i + sampler.getNext1D()) - height) / height;
        Ray ray;
        float pdf;
        if (camera.sampleRay(Vec2(u, v), ray, pdf)) {
          const Vec3 radiance = integrator.integrate(ray, scene, sampler) / pdf;
          image.addPixel(i, j, radiance);
        } else {
          image.setPixel(i, j, Vec3(0));
        }
      }
    }
  }

  // take average
  image.divide(n_samples);

  image.gammaCorrection(2.2f);
  image.writePPM("output.ppm");
}