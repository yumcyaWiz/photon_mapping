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
  const int n_photons = 1000000;
  const int max_depth = 100;
  const Vec3 camPos(2.78, 2.73, -9);
  const Vec3 lookAt(2.78, 2.73, 2.796);

  Image image(width, height);

  const Camera camera(camPos, normalize(lookAt - camPos), 0.25 * PI);

  // cornellbox scene
  const auto white = std::make_shared<Lambert>(Vec3(0.8));
  const auto red = std::make_shared<Lambert>(Vec3(0.8, 0.05, 0.05));
  const auto green = std::make_shared<Lambert>(Vec3(0.05, 0.8, 0.05));

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
  scene.addPrimitive(Primitive(floor, white));
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
  PhotonMapping integrator(n_photons, 1, max_depth);
  UniformSampler sampler;
  integrator.build(scene, sampler);

  // visualize photon map
  spdlog::info("[main] visualizing photon map");

  const PhotonMap* photon_map = integrator.getPhotonMapPtr();

#pragma omp parallel for collapse(2)
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      const float u = (2.0f * j - width) / height;
      const float v = (2.0f * i - height) / height;
      Ray ray;
      float pdf;
      if (camera.sampleRay(Vec2(u, v), ray, pdf)) {
        IntersectInfo info;
        if (scene.intersect(ray, info)) {
          // query photon map
          float r2;
          const int photon_idx =
              photon_map->queryNearestPhoton(info.surfaceInfo.position, r2);

          // if distance to the photon is small enough, write photon's
          // throughput to the image
          if (r2 < 0.001f) {
            const Photon& photon = photon_map->getIthPhoton(photon_idx);
            image.setPixel(i, j, photon.throughput);
          }
        } else {
          image.setPixel(i, j, Vec3(0));
        }
      } else {
        image.setPixel(i, j, Vec3(0));
      }
    }
  }

  image.writePPM("output.ppm");

  return 0;
}