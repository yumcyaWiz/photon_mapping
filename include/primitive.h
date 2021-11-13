#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H
#include <cmath>
#include <memory>

#include "light.h"
#include "material.h"
#include "shape.h"

class Primitive {
 private:
  const std::shared_ptr<Shape> shape;
  const std::shared_ptr<BxDF> bxdf;
  const std::shared_ptr<AreaLight> areaLight;

 public:
  Primitive(const std::shared_ptr<Shape>& shape,
            const std::shared_ptr<BxDF>& bxdf,
            const std::shared_ptr<AreaLight>& areaLight = nullptr)
      : shape(shape), bxdf(bxdf), areaLight(areaLight) {}

  bool hasAreaLight() const { return areaLight != nullptr; }

  // NOTE: for populating scene's light array
  std::shared_ptr<AreaLight> getAreaLightPtr() const { return areaLight; }

  Vec3 Le(const SurfaceInfo& surfInfo, const Vec3& dir) const {
    return areaLight->Le(surfInfo, dir);
  }

  bool intersect(const Ray& ray, IntersectInfo& info) const {
    if (shape->intersect(ray, info)) {
      info.hitPrimitive = this;
      return true;
    }
    return false;
  }

  BxDFType getBxDFType() const { return bxdf->getType(); }

  Vec3 sampleBRDF(const Vec3& wo, const SurfaceInfo& surfInfo, Sampler& sampler,
                  Vec3& wi, float& pdf) const {
    // world to local transform
    const Vec3 wo_l =
        worldToLocal(wo, surfInfo.dpdu, surfInfo.normal, surfInfo.dpdv);

    // sample direction in tangent space
    Vec3 wi_l;
    const Vec3 f = bxdf->sampleDirection(wo_l, sampler, wi_l, pdf);

    // local to world transform
    wi = localToWorld(wi_l, surfInfo.dpdu, surfInfo.normal, surfInfo.dpdv);

    return f;
  }
};

#endif