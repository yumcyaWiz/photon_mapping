#include <omp.h>

#include "camera.h"
#include "image.h"
#include "integrator.h"
#include "photon_map.h"
#include "scene.h"

int main() {
  const int width = 512;
  const int height = 512;
  const int n_photons = 1000000;
  const int max_depth = 100;

  Image image(width, height);
  Camera camera(Vec3f(0, 1, 7), Vec3f(0, 0, -1), 0.25 * PI);

  Scene scene;
  scene.loadModel("CornellBox-Original.obj");
  scene.build();

  // photon tracing and build photon map
  PhotonMapping integrator(n_photons, 1, 0, 0, 0, false, max_depth);
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
      if (camera.sampleRay(Vec2f(u, v), ray, pdf)) {
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
          image.setPixel(i, j, Vec3f(0));
        }
      } else {
        image.setPixel(i, j, Vec3f(0));
      }
    }
  }

  image.writePPM("output.ppm");

  return 0;
}