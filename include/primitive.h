#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H
#include <cmath>
#include <memory>

#include "light.h"
#include "material.h"
#include "triangle.h"

// primitive provides an abstraction layer of the object's shape(triangle),
// material, area light
class Primitive {
 private:
  const Triangle* triangle;
  const std::shared_ptr<BxDF> bxdf;
  const std::shared_ptr<Light> areaLight;

 public:
  Primitive(const Triangle* triangle, const std::shared_ptr<BxDF>& bxdf,
            const std::shared_ptr<Light>& areaLight = nullptr)
      : triangle(triangle), bxdf(bxdf), areaLight(areaLight) {}

  bool hasAreaLight() const { return areaLight != nullptr; }

  // return emission
  Vec3f Le(const SurfaceInfo& surfInfo, const Vec3f& dir) const {
    return areaLight->Le(surfInfo, dir);
  }

  BxDFType getBxDFType() const { return bxdf->getType(); }

  Vec3f evaluateBxDF(const Vec3f& wo, const Vec3f& wi,
                     const SurfaceInfo& surfInfo,
                     const TransportDirection& mode) const {
    // world to local transform
    const Vec3f wo_l =
        worldToLocal(wo, surfInfo.dpdu, surfInfo.shadingNormal, surfInfo.dpdv);
    const Vec3f wi_l =
        worldToLocal(wi, surfInfo.dpdu, surfInfo.shadingNormal, surfInfo.dpdv);

    return bxdf->evaluate(wo_l, wi_l, mode);
  }

  // sample direction by BxDF
  // its pdf is propotional to the shape od BxDF
  Vec3f sampleBxDF(const Vec3f& wo, const SurfaceInfo& surfInfo,
                   const TransportDirection& mode, Sampler& sampler, Vec3f& wi,
                   float& pdf) const {
    // world to local transform
    const Vec3f wo_l =
        worldToLocal(wo, surfInfo.dpdu, surfInfo.shadingNormal, surfInfo.dpdv);

    // sample direction in tangent space
    Vec3f wi_l;
    const Vec3f f = bxdf->sampleDirection(wo_l, mode, sampler, wi_l, pdf);

    // local to world transform
    wi = localToWorld(wi_l, surfInfo.dpdu, surfInfo.shadingNormal,
                      surfInfo.dpdv);

    return f;
  }

  // get all samplable direction
  std::vector<DirectionPair> sampleAllBxDF(
      const Vec3f& wo, const SurfaceInfo& surfInfo,
      const TransportDirection& mode) const {
    // world to local transform
    const Vec3f wo_l =
        worldToLocal(wo, surfInfo.dpdu, surfInfo.shadingNormal, surfInfo.dpdv);

    // sample all direction in tangent space
    std::vector<DirectionPair> dir_pairs = bxdf->sampleAllDirection(wo_l, mode);

    // local to world transform
    for (auto& dp : dir_pairs) {
      dp.first = localToWorld(dp.first, surfInfo.dpdu, surfInfo.shadingNormal,
                              surfInfo.dpdv);
    }

    return dir_pairs;
  }
};

#endif