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
  const std::shared_ptr<Material> material;
  const std::shared_ptr<AreaLight> areaLight;

 public:
  Primitive(const std::shared_ptr<Shape>& shape,
            const std::shared_ptr<Material>& material,
            const std::shared_ptr<AreaLight>& areaLight = nullptr)
      : shape(shape), material(material), areaLight(areaLight) {}

  bool hasAreaLight() const { return areaLight != nullptr; }

  std::shared_ptr<AreaLight> getAreaLightPtr() const { return areaLight; }

  bool intersect(const Ray& ray, IntersectInfo& info) const {
    if (shape->intersect(ray, info)) {
      info.hitPrimitive = this;
      return true;
    }
    return false;
  }

  Vec3 sampleBRDF(const Vec3& wo, const SurfaceInfo& surfInfo, Sampler& sampler,
                  Vec3& wi, MaterialType& type, float& pdf) const {
    // world to local transform
    const Vec3 wo_l =
        worldToLocal(wo, surfInfo.dpdu, surfInfo.normal, surfInfo.dpdv);

    // sample direction in tangent space
    Vec3 wi_l;
    const Vec3 brdf = material->sampleBRDF(wo_l, sampler, wi_l, type, pdf);

    // local to world transform
    wi = localToWorld(wi_l, surfInfo.dpdu, surfInfo.normal, surfInfo.dpdv);

    return brdf;
  }
};

#endif