#ifndef _MATERIAL_H
#define _MATERIAL_H
#include <memory>

#include "core.h"
#include "sampler.h"

enum class BxDFType { DIFFUSE, SPECULAR };

class BxDF {
 private:
  BxDFType type;

 public:
  BxDF(const BxDFType& type) : type(type) {}

  static float cosTheta(const Vec3& v) { return v[1]; }

  BxDFType getType() const { return type; }

  virtual Vec3 evaluate(const Vec3& wo, const Vec3& wi) const = 0;
  virtual Vec3 sampleDirection(const Vec3& wo, Sampler& sampler, Vec3& wi,
                               float& pdf) const = 0;
};

class Lambert : public BxDF {
 private:
  Vec3 rho;

 public:
  Lambert(const Vec3& rho) : BxDF(BxDFType::DIFFUSE), rho(rho) {}

  Vec3 evaluate(const Vec3& wo, const Vec3& wi) const {
    // when wo, wi is under the surface, return 0
    const float cosThetaO = BxDF::cosTheta(wo);
    const float cosThetaI = BxDF::cosTheta(wi);
    if (cosThetaO < 0 || cosThetaI < 0) return Vec3(0);

    return rho / PI;
  }

  Vec3 sampleDirection(const Vec3& wo, Sampler& sampler, Vec3& wi,
                       float& pdf) const override {
    // cosine weighted hemisphere sampling
    wi = sampleCosineHemisphere(sampler.getNext2D(), pdf);

    return evaluate(wo, wi);
  }
};

#endif