#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H
#include "core.h"
#include "photon_map.h"
#include "scene.h"

class Integrator {
 public:
  // do preliminary jobs before calling integrate
  virtual void build(const Scene& scene, Sampler& sampler) const = 0;

  // compute radiance coming from given ray
  virtual Vec3 integrate(const Ray& ray, const Scene& scene,
                         Sampler& sampler) const = 0;
};

// implementation of photon mapping
class PM : public Integrator {
 private:
  const int nPhotons;
  const int nDensityEstimation;
  static constexpr int maxDepth = 100;
  PhotonMap photonMap;

 public:
  PM(int nPhotons, int nDensityEstimation)
      : nPhotons(nPhotons), nDensityEstimation(nDensityEstimation) {}

  void build(const Scene& scene, Sampler& sampler) const override {
    // TODO: trace photons, build photon map
    std::vector<Photon> photons(nPhotons);

    for (int i = 0; i < nPhotons; ++i) {
      // sample light
      float light_choose_pdf;
      const std::shared_ptr<Light> light =
          scene.sampleLight(sampler, light_choose_pdf);

      // sample point on light
      float light_pos_pdf;
      const SurfaceInfo lightSurf = light->samplePoint(sampler, light_pos_pdf);

      // sample direction on light
      float light_dir_pdf;
      const Vec3 lightDir =
          light->sampleDirection(lightSurf, sampler, light_dir_pdf);

      // spawn ray
      Ray ray(lightSurf.position, light_dir_pdf);

      // trace photons until hitting diffuse surface
      for (int k = 0; k < maxDepth; ++k) {
      }
    }
  }

  Vec3 integrate(const Ray& ray, const Scene& scene,
                 Sampler& sampler) const override {
    // TODO: trace from camera until hitting diffuse surface, compute radiance
    return Vec3(0);
  }
};

#endif