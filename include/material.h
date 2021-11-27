#ifndef _MATERIAL_H
#define _MATERIAL_H
#include <memory>

#include "core.h"
#include "sampler.h"

enum class BxDFType { DIFFUSE, SPECULAR };

enum class TransportDirection { FROM_LIGHT, FROM_CAMERA };

using DirectionPair = std::pair<Vec3f, Vec3f>;

// represent BRDF or BTDF
// direction vectors are in tangent space(x: tangent, y: normal, z: bitangent)
class BxDF {
 private:
  BxDFType type;

 public:
  BxDF(const BxDFType& type) : type(type) {}

  static float cosTheta(const Vec3f& v) { return v[1]; }
  static float absCosTheta(const Vec3f& v) { return std::abs(cosTheta(v)); }

  // compute reflection direction
  static Vec3f reflect(const Vec3f& v, const Vec3f& n) {
    return -v + 2.0f * dot(v, n) * n;
  }

  // compute refracted direction
  static bool refract(const Vec3f& v, const Vec3f& n, float iorI, float iorT,
                      Vec3f& t) {
    const Vec3f t_h = -iorI / iorT * (v - dot(v, n) * n);
    // total reflection
    if (length2(t_h) > 1.0f) {
      return false;
    }
    const Vec3f t_p = -std::sqrt(std::max(1.0f - length2(t_h), 0.0f)) * n;
    t = t_h + t_p;
    return true;
  }

  // schlick approximation of fresnel reflectance
  static float fresnel(float cosThetaI, float iorI, float iorT) {
    const float f0 =
        (iorI - iorT) * (iorI - iorT) / ((iorI + iorT) * (iorI + iorT));
    const auto pow5 = [](float x) { return x * x * x * x * x; };
    return f0 + (1.0f - f0) * pow5(std::max(1.0f - std::abs(cosThetaI), 0.0f));
  }

  // get BxDF type
  BxDFType getType() const { return type; }

  // evaluate BxDF
  virtual Vec3f evaluate(const Vec3f& wo, const Vec3f& wi,
                         const TransportDirection& transport_dir) const = 0;

  // sample direction by BxDF.
  // its pdf is propotional to the shape of BxDF
  virtual Vec3f sampleDirection(const Vec3f& wo,
                                const TransportDirection& transport_dir,
                                Sampler& sampler, Vec3f& wi,
                                float& pdf) const = 0;

  // get all samplable direction
  // NOTE: for specular only
  // NOTE: used for drawing fresnel reflection nicely at low number of samples
  virtual std::vector<DirectionPair> sampleAllDirection(
      const Vec3f& wo, const TransportDirection& transport_dir) const = 0;
};

class Lambert : public BxDF {
 private:
  Vec3f rho;

 public:
  Lambert(const Vec3f& rho) : BxDF(BxDFType::DIFFUSE), rho(rho) {}

  Vec3f evaluate(const Vec3f& wo, const Vec3f& wi,
                 const TransportDirection& transport_dir) const override {
    // when wo, wi is under the surface, return 0
    const float cosThetaO = cosTheta(wo);
    const float cosThetaI = cosTheta(wi);
    if (cosThetaO < 0 || cosThetaI < 0) return Vec3f(0);

    return rho / PI;
  }

  Vec3f sampleDirection(const Vec3f& wo,
                        const TransportDirection& transport_dir,
                        Sampler& sampler, Vec3f& wi,
                        float& pdf) const override {
    // cosine weighted hemisphere sampling
    wi = sampleCosineHemisphere(sampler.getNext2D(), pdf);

    return evaluate(wo, wi, transport_dir);
  }

  std::vector<DirectionPair> sampleAllDirection(
      const Vec3f& wo, const TransportDirection& transport_dir) const override {
    std::vector<DirectionPair> ret;
    return ret;
  }
};

class Mirror : public BxDF {
 private:
  Vec3f rho;

 public:
  Mirror(const Vec3f& rho) : BxDF(BxDFType::SPECULAR), rho(rho) {}

  // NOTE: delta function
  Vec3f evaluate(const Vec3f& wo, const Vec3f& wi,
                 const TransportDirection& transport_dir) const override {
    return Vec3f(0);
  }

  Vec3f sampleDirection(const Vec3f& wo,
                        const TransportDirection& transport_dir,
                        Sampler& sampler, Vec3f& wi,
                        float& pdf) const override {
    wi = reflect(wo, Vec3f(0, 1, 0));
    pdf = 1.0f;

    return rho / absCosTheta(wi);
  }

  std::vector<DirectionPair> sampleAllDirection(
      const Vec3f& wo, const TransportDirection& transport_dir) const override {
    std::vector<DirectionPair> ret;
    const Vec3f wi = reflect(wo, Vec3f(0, 1, 0));
    ret.emplace_back(wi, rho / absCosTheta(wi));
    return ret;
  }
};

// NOTE: due to the asymmetry of BSDF, we need to use a different scaling factor
// for photon tracing
// https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation#Non-symmetricScattering
class Glass : public BxDF {
 private:
  Vec3f rho;
  float ior;

 public:
  Glass(const Vec3f& rho, float ior)
      : BxDF(BxDFType::SPECULAR), rho(rho), ior(ior) {}

  // NOTE: delta function
  Vec3f evaluate(const Vec3f& wo, const Vec3f& wi,
                 const TransportDirection& transport_dir) const override {
    return Vec3f(0);
  }

  Vec3f sampleDirection(const Vec3f& wo,
                        const TransportDirection& transport_dir,
                        Sampler& sampler, Vec3f& wi,
                        float& pdf) const override {
    // set appropriate ior, normal
    float iorO, iorI;
    Vec3f n;
    if (wo[1] > 0) {
      iorO = 1.0f;
      iorI = ior;
      n = Vec3f(0, 1, 0);
    } else {
      iorO = ior;
      iorI = 1.0f;
      n = Vec3f(0, -1, 0);
    }

    // fresnel reflectance
    const float fr = fresnel(dot(wo, n), iorO, iorI);

    // reflection
    if (sampler.getNext1D() < fr) {
      wi = reflect(wo, n);
      pdf = 1.0f;
      return rho / absCosTheta(wi);
    }
    // refraction
    else {
      Vec3f tr;
      if (refract(wo, n, iorO, iorI, tr)) {
        wi = tr;
        pdf = 1.0f;

        float scalling = 1.0f;
        if (transport_dir == TransportDirection::FROM_CAMERA) {
          scalling = (iorO * iorO) / (iorI * iorI);
        }

        return scalling * rho / absCosTheta(wi);
      }
      // total reflection
      else {
        wi = reflect(wo, n);
        pdf = 1.0f;
        return rho / absCosTheta(wi);
      }
    }
  }

  std::vector<DirectionPair> sampleAllDirection(
      const Vec3f& wo, const TransportDirection& transport_dir) const override {
    std::vector<DirectionPair> ret;

    // set appropriate ior, normal
    float iorO, iorI;
    Vec3f n;
    if (wo[1] > 0) {
      iorO = 1.0f;
      iorI = ior;
      n = Vec3f(0, 1, 0);
    } else {
      iorO = ior;
      iorI = 1.0f;
      n = Vec3f(0, -1, 0);
    }

    // fresnel reflectance
    const float fr = fresnel(dot(wo, n), iorO, iorI);

    // reflection
    const Vec3f wr = reflect(wo, n);
    ret.emplace_back(wr, fr * rho / absCosTheta(wr));

    // refraction
    Vec3f tr;
    if (refract(wo, n, iorO, iorI, tr)) {
      float scalling = 1.0f;
      if (transport_dir == TransportDirection::FROM_CAMERA) {
        scalling = (iorO * iorO) / (iorI * iorI);
      }

      ret.emplace_back(tr, (1.0f - fr) * scalling * rho / absCosTheta(tr));
    } else {
      ret[0].second = rho / absCosTheta(wr);
    }

    return ret;
  }
};

#endif