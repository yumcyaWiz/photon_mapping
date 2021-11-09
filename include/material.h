#ifndef _MATERIAL_H
#define _MATERIAL_H
#include "core.h"
#include "sampler.h"

class Material {
  virtual Vec3 BRDF(const Vec3& wo, const Vec3& wi) const = 0;
  virtual Vec3 sampleBRDF(const Vec3& wo, Sampler& sampler, Vec3& wi,
                          float& pdf) const;
};

class Diffuse : public Material {
 private:
  Vec3 rho;

  Diffuse(const Vec3& rho) : rho(rho) {}

  Vec3 BRDF(const Vec3& wo, const Vec3& wi) const { return rho / PI; }

  Vec3 sampleBRDF(const Vec3& wo, Sampler& sampler, Vec3& wi,
                  float& pdf) const override {
    // cosine weighted hemisphere sampling
    wi = sampleCosineHemisphere(sampler.getNext2D(), pdf);
    return BRDF(wo, wi);
  }
};

class Mirror : public Material {};

class Glass : public Material {};

#endif