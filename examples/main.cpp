#include "camera.h"
#include "image.h"
#include "integrator.h"
#include "photon_map.h"
#include "scene.h"

int main() {
  const int width = 512;
  const int height = 512;
  const int n_samples = 100;
  const int n_photons = 100000;
  const int n_estimation_global = 100;
  const float n_photons_caustics_multiplier = 100;
  const int n_estimation_caustics = 100;
  const bool final_gathering = true;
  const int strict_calc_depth = 3;
  const int max_depth = 100;

  Image image(width, height);
  Camera camera(Vec3f(0, 1, 6), Vec3f(0, 0, -1), 0.25 * PI);

  Scene scene;
  scene.loadModel("CornellBox-Water.obj");
  scene.build();

  // photon tracing and build photon map
  PhotonMapping integrator(n_photons, n_estimation_global,
                           n_photons_caustics_multiplier, n_estimation_caustics,
                           strict_calc_depth, final_gathering, max_depth);
  UniformSampler sampler;
  integrator.build(scene, sampler);

  spdlog::info("[main] tracing rays from camera");
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
        if (camera.sampleRay(Vec2f(u, v), ray, pdf)) {
          const Vec3f radiance =
              integrator.integrate(ray, scene, sampler) / pdf;

          if (std::isnan(radiance[0]) || std::isnan(radiance[1]) ||
              std::isnan(radiance[2])) {
            spdlog::error("radiance is NaN");
            continue;
          } else if (radiance[0] < 0 || radiance[1] < 0 || radiance[2] < 0) {
            spdlog::error("radiance is minus");
            continue;
          }

          image.addPixel(i, j, radiance);
        } else {
          image.setPixel(i, j, Vec3f(0));
        }
      }
    }
  }

  // take average
  image.divide(n_samples);

  image.gammaCorrection(2.2f);
  image.writePPM("output.ppm");
}