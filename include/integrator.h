#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H
#include <optional>

#include "core.h"
#include "photon_map.h"
#include "scene.h"

class Integrator {
 public:
  // do preliminary jobs before calling integrate
  virtual void build(const Scene& scene, Sampler& sampler) = 0;

  // compute radiance coming from given ray
  virtual Vec3 integrate(const Ray& ray, const Scene& scene,
                         Sampler& sampler) const = 0;
};

// implementation of photon mapping
class PhotonMapping : public Integrator {
 private:
  const int nPhotons;
  const int nDensityEstimation;
  static constexpr int maxDepth = 100;
  PhotonMap photon_map;

 public:
  PhotonMapping(int nPhotons, int nDensityEstimation)
      : nPhotons(nPhotons), nDensityEstimation(nDensityEstimation) {}

  const PhotonMap* getPhotonMapPtr() const { return &photon_map; }

  void build(const Scene& scene, Sampler& sampler) override {
    // NOTE: to parallelize photon tracing, prepare array of std::option<Photon>
    std::vector<std::optional<Photon>> photons(nPhotons);

    // photon tracing
    spdlog::info("[PhotonMapping] tracing photons");
    for (int i = 0; i < nPhotons; ++i) {
      // sample light
      float light_choose_pdf;
      const std::shared_ptr<Light> light =
          scene.sampleLight(sampler, light_choose_pdf);

      // sample point on light
      float light_pos_pdf;
      const SurfaceInfo light_surf = light->samplePoint(sampler, light_pos_pdf);

      // sample direction on light
      float light_dir_pdf;
      const Vec3 light_dir =
          light->sampleDirection(light_surf, sampler, light_dir_pdf);

      spdlog::info("{}, {}, {}", light_surf.normal[0], light_surf.normal[1],
                   light_surf.normal[2]);

      // spawn ray
      Ray ray(light_surf.position, light_dir_pdf);
      Vec3 throughput = light->Le(light_surf, light_dir) /
                        (light_choose_pdf * light_pos_pdf) *
                        dot(light_dir, light_surf.normal);

      // trace photons
      // whener hitting diffuse surface, add photon to the photon array
      // recursively tracing photon with russian roulette
      for (int k = 0; k < 1; ++k) {
        IntersectInfo info;
        if (scene.intersect(ray, info)) {
          // if hitting diffuse surface, add photon to the photon array
          BxDFType bxdf_type = info.hitPrimitive->getBxDFType();
          if (bxdf_type == BxDFType::DIFFUSE) {
            photons[i] =
                Photon(throughput, info.surfaceInfo.position, -ray.direction);
          }

          // russian roulette
          if (k > 0) {
            const float russian_roulette_prob = std::min(
                std::max(throughput[0], std::max(throughput[1], throughput[2])),
                1.0f);
            if (sampler.getNext1D() >= russian_roulette_prob) {
              break;
            }
            throughput /= russian_roulette_prob;
          }

          // sample direction by BxDF
          Vec3 dir;
          float pdf_dir;
          Vec3 f = info.hitPrimitive->sampleBRDF(
              -ray.direction, info.surfaceInfo, sampler, dir, pdf_dir);

          // update throughput and ray
          throughput *=
              f * std::abs(dot(dir, info.surfaceInfo.normal)) / pdf_dir;
          ray = Ray(info.surfaceInfo.position, dir);
        } else {
          // photon goes out to the sky
          break;
        }
      }
    }

    // add photon to the photon map
    for (const auto& p : photons) {
      if (p) {
        photon_map.addPhoton(p.value());
      }
    }

    // build photon map
    spdlog::info("[PhotonMapping] building photon map");
    photon_map.build();
  }

  Vec3 integrate(const Ray& ray, const Scene& scene,
                 Sampler& sampler) const override {
    // TODO: trace from camera until hitting diffuse surface, compute radiance
    return Vec3(0);
  }
};

#endif