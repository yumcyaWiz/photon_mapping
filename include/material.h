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

  BxDFType getType() const { return type; }

  virtual Vec3 evaluate(const Vec3& wo, const Vec3& wi) const = 0;
  virtual Vec3 sampleDirection(const Vec3& wo, Sampler& sampler, Vec3& wi,
                               float& pdf) const;
};

class Lambert : public BxDF {
 private:
  Vec3 rho;

 public:
  Lambert(const Vec3& rho) : BxDF(BxDFType::DIFFUSE), rho(rho) {}

  Vec3 evaluate(const Vec3& wo, const Vec3& wi) const { return rho / PI; }

  Vec3 sampleDirection(const Vec3& wo, Sampler& sampler, Vec3& wi,
                       float& pdf) const override {
    // cosine weighted hemisphere sampling
    wi = sampleCosineHemisphere(sampler.getNext2D(), pdf);

    return evaluate(wo, wi);
  }
};

#endif