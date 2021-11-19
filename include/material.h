#ifndef _MATERIAL_H
#define _MATERIAL_H
#include <memory>

#include "core.h"
#include "sampler.h"

enum class BxDFType { DIFFUSE, SPECULAR };

using DirectionPair = std::pair<Vec3, Vec3>;

class BxDF {
 private:
  BxDFType type;

 public:
  BxDF(const BxDFType& type) : type(type) {}

  static float cosTheta(const Vec3& v) { return v[1]; }
  static float absCosTheta(const Vec3& v) { return std::abs(cosTheta(v)); }

  static Vec3 reflect(const Vec3& v, const Vec3& n) {
    return -v + 2.0f * dot(v, n) * n;
  }
  static bool refract(const Vec3& v, const Vec3& n, float iorI, float iorT,
                      Vec3& t) {
    const Vec3 t_h = -iorI / iorT * (v - dot(v, n) * n);
    // total reflection
    if (length2(t_h) > 1.0f) {
      return false;
    }
    const Vec3 t_p = -std::sqrt(std::max(1.0f - length2(t_h), 0.0f)) * n;
    t = t_h + t_p;
    return true;
  }

  // schlick approximation
  static float fresnel(float cosThetaI, float iorI, float iorT) {
    const float f0 =
        (iorI - iorT) * (iorI - iorT) / ((iorI + iorT) * (iorI + iorT));
    const auto pow5 = [](float x) { return x * x * x * x * x; };
    return f0 + (1.0f - f0) * pow5(std::max(1.0f - std::abs(cosThetaI), 0.0f));
  }

  // get BxDF type
  BxDFType getType() const { return type; }

  // evaluate BxDF
  virtual Vec3 evaluate(const Vec3& wo, const Vec3& wi) const = 0;

  // sample direction
  // pdf is set to be propotional to BxDF
  virtual Vec3 sampleDirection(const Vec3& wo, Sampler& sampler, Vec3& wi,
                               float& pdf) const = 0;

  // sample all samplable direction
  // NOTE: for specular only
  // NOTE: used for drawing fresnel reflection at low number of samples
  virtual std::vector<DirectionPair> sampleAllDirection(
      const Vec3& wo) const = 0;
};

class Lambert : public BxDF {
 private:
  Vec3 rho;

 public:
  Lambert(const Vec3& rho) : BxDF(BxDFType::DIFFUSE), rho(rho) {}

  Vec3 evaluate(const Vec3& wo, const Vec3& wi) const override {
    // when wo, wi is under the surface, return 0
    const float cosThetaO = cosTheta(wo);
    const float cosThetaI = cosTheta(wi);
    if (cosThetaO < 0 || cosThetaI < 0) return Vec3(0);

    return rho / PI;
  }

  Vec3 sampleDirection(const Vec3& wo, Sampler& sampler, Vec3& wi,
                       float& pdf) const override {
    // cosine weighted hemisphere sampling
    wi = sampleCosineHemisphere(sampler.getNext2D(), pdf);

    return evaluate(wo, wi);
  }

  std::vector<DirectionPair> sampleAllDirection(const Vec3& wo) const override {
    std::vector<DirectionPair> ret;
    return ret;
  }
};

class Mirror : public BxDF {
 private:
  Vec3 rho;

 public:
  Mirror(const Vec3& rho) : BxDF(BxDFType::SPECULAR), rho(rho) {}

  // NOTE: delta function
  Vec3 evaluate(const Vec3& wo, const Vec3& wi) const override {
    return Vec3(0);
  }

  Vec3 sampleDirection(const Vec3& wo, Sampler& sampler, Vec3& wi,
                       float& pdf) const override {
    wi = reflect(wo, Vec3(0, 1, 0));
    pdf = 1.0f;

    return rho / absCosTheta(wi);
  }

  std::vector<DirectionPair> sampleAllDirection(const Vec3& wo) const override {
    std::vector<DirectionPair> ret;
    const Vec3 wi = reflect(wo, Vec3(0, 1, 0));
    ret.emplace_back(wi, rho / absCosTheta(wi));
    return ret;
  }
};

class Glass : public BxDF {
 private:
  Vec3 rho;
  float ior;

 public:
  Glass(const Vec3& rho, float ior)
      : BxDF(BxDFType::SPECULAR), rho(rho), ior(ior) {}

  // NOTE: delta function
  Vec3 evaluate(const Vec3& wo, const Vec3& wi) const override {
    return Vec3(0);
  }

  Vec3 sampleDirection(const Vec3& wo, Sampler& sampler, Vec3& wi,
                       float& pdf) const override {
    // set appropriate ior, normal
    float iorI, iorT;
    Vec3 n;
    if (wo[1] > 0) {
      iorI = 1.0f;
      iorT = ior;
      n = Vec3(0, 1, 0);
    } else {
      iorI = ior;
      iorT = 1.0f;
      n = Vec3(0, -1, 0);
    }

    // fresnel reflectance
    const float fr = fresnel(dot(wo, n), iorI, iorT);

    // reflection
    if (sampler.getNext1D() < fr) {
      wi = reflect(wo, n);
      pdf = 1.0f;
      return rho / absCosTheta(wi);
    }
    // refraction
    else {
      Vec3 tr;
      if (refract(wo, n, iorI, iorT, tr)) {
        wi = tr;
        pdf = 1.0f;
        return iorT * iorT / (iorI * iorI) * rho / absCosTheta(wi);
      }
      // total reflection
      else {
        wi = reflect(wo, n);
        pdf = 1.0f;
        return rho / absCosTheta(wi);
      }
    }
  }

  std::vector<DirectionPair> sampleAllDirection(const Vec3& wo) const override {
    std::vector<DirectionPair> ret;

    // set appropriate ior, normal
    float iorI, iorT;
    Vec3 n;
    if (wo[1] > 0) {
      iorI = 1.0f;
      iorT = ior;
      n = Vec3(0, 1, 0);
    } else {
      iorI = ior;
      iorT = 1.0f;
      n = Vec3(0, -1, 0);
    }

    // fresnel reflectance
    const float fr = fresnel(dot(wo, n), iorI, iorT);

    // reflection
    const Vec3 wr = reflect(wo, n);
    ret.emplace_back(wr, fr * rho / absCosTheta(wr));

    // refraction
    Vec3 tr;
    if (refract(wo, n, iorI, iorT, tr)) {
      ret.emplace_back(tr, std::max((1.0f - fr), 0.0f) * (iorT * iorT) /
                               (iorI * iorI) * rho / absCosTheta(tr));
    } else {
      ret[0].second = rho / absCosTheta(wr);
    }

    return ret;
  }
};

#endif