#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H
#include "core.h"
#include "photon_map.h"
#include "scene.h"

class Integrator {
 public:
  virtual void build(const Scene& scene, Sampler& sampler) const = 0;
  virtual Vec3 integrate(const Ray& ray, const Scene& scene,
                         Sampler& sampler) const = 0;
};

// implementation of photon mapping
class PM : public Integrator {
 private:
  const int nPhotons;
  const int nDensityEstimation;
  PhotonMap photonMap;

 public:
  PM(int nPhotons, int nDensityEstimation)
      : nPhotons(nPhotons), nDensityEstimation(nDensityEstimation) {}

  void build(const Scene& scene, Sampler& sampler) const override {
    // TODO: trace photons, build photon map
  }

  Vec3 integrate(const Ray& ray, const Scene& scene,
                 Sampler& sampler) const override {
    // TODO: trace until hitting diffuse surface, compute radiance
    return Vec3(0);
  }
};

#endif