#ifndef _LIGHT_H
#define _LIGHT_H
#include <memory>

#include "core.h"
#include "primitive.h"
#include "sampler.h"
#include "shape.h"

class Light {
 public:
  virtual Vec3 samplePoint(Sampler& sampler, Vec3& pos, Vec3& dir,
                           float& pdf) const = 0;
};

class AreaLight : public Light {
 private:
  const Vec3 le;
  const std::shared_ptr<Shape> shape;

 public:
  AreaLight(const Vec3& le, const std::shared_ptr<Shape>& shape)
      : le(le), shape(shape) {}

  Vec3 samplePoint(Sampler& sampler, Vec3& pos, Vec3& dir,
                   float& pdf) const override {
    // sample point on shape
    float pdf_point;
    SurfaceInfo surfInfo = shape->samplePoint(sampler, pdf_point);
    pos = surfInfo.position;

    // sample direction by cosine weighted hemisphere sampling
    float pdf_dir;
    dir = sampleCosineHemisphere(sampler.getNext2D(), pdf_dir);

    return le;
  }
};

#endif