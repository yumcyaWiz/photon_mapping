#ifndef _LIGHT_H
#define _LIGHT_H
#include <memory>

#include "core.h"
#include "sampler.h"
#include "shape.h"

class Light {
 public:
  virtual Vec3 Le(const SurfaceInfo& info, const Vec3& dir) const = 0;
  virtual SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const = 0;
  virtual Vec3 sampleDirection(const SurfaceInfo& surfInfo, Sampler& sampler,
                               float& pdf) const = 0;
};

class AreaLight : public Light {
 private:
  const Vec3 le;
  const std::shared_ptr<Shape> shape;

 public:
  AreaLight(const Vec3& le, const std::shared_ptr<Shape>& shape)
      : le(le), shape(shape) {}

  Vec3 Le(const SurfaceInfo& info, const Vec3& dir) const override {
    return le;
  }

  SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const override {
    return shape->samplePoint(sampler, pdf);
  }

  Vec3 sampleDirection(const SurfaceInfo& surfInfo, Sampler& sampler,
                       float& pdf) const override {
    const Vec3 dir = sampleCosineHemisphere(sampler.getNext2D(), pdf);

    // transform direction from local to world
    return localToWorld(dir, surfInfo.dpdu, surfInfo.normal, surfInfo.dpdv);
  }
};

#endif