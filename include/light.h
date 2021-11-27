#ifndef _LIGHT_H
#define _LIGHT_H

#include "core.h"
#include "sampler.h"
#include "triangle.h"

// light interface
class Light {
 public:
  virtual Vec3f Le(const SurfaceInfo& info, const Vec3f& dir) const = 0;
  virtual SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const = 0;
  virtual Vec3f sampleDirection(const SurfaceInfo& surfInfo, Sampler& sampler,
                                float& pdf) const = 0;
};

class AreaLight : public Light {
 private:
  const Vec3f le;  // emission
  const Triangle* triangle;

 public:
  AreaLight(const Vec3f& le, const Triangle* triangle)
      : le(le), triangle(triangle) {}

  // return emission
  Vec3f Le(const SurfaceInfo& info, const Vec3f& dir) const override {
    return le;
  }

  // sample point on the light
  SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const override {
    return triangle->samplePoint(sampler, pdf);
  }

  // sample direction from the light
  Vec3f sampleDirection(const SurfaceInfo& surfInfo, Sampler& sampler,
                        float& pdf) const override {
    const Vec3f dir = sampleCosineHemisphere(sampler.getNext2D(), pdf);

    // transform direction from local to world
    return localToWorld(dir, surfInfo.dpdu, surfInfo.shadingNormal,
                        surfInfo.dpdv);
  }
};

#endif