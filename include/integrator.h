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

// implementation of path tracing
// for reference purpose
class PathTracing : public Integrator {
 private:
  const int maxDepth;

 public:
  PathTracing(int maxDepth = 100) : maxDepth(maxDepth) {}

  void build(const Scene& scene, Sampler& sampler) override {}

  Vec3 integrate(const Ray& ray_in, const Scene& scene,
                 Sampler& sampler) const override {
    Ray ray = ray_in;
    Vec3 throughput(1, 1, 1);

    for (int k = 0; k < maxDepth; ++k) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
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

        // Le
        if (info.hitPrimitive->hasAreaLight()) {
          return throughput *
                 info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
        }

        // sample direction by BxDF
        Vec3 dir;
        float pdf_dir;
        Vec3 f = info.hitPrimitive->sampleBxDF(-ray.direction, info.surfaceInfo,
                                               sampler, dir, pdf_dir);

        // update throughput and ray
        throughput *= f * std::abs(dot(dir, info.surfaceInfo.normal)) / pdf_dir;
        ray = Ray(info.surfaceInfo.position, dir);
      } else {
        break;
      }
    }

    return Vec3(0);
  }
};

// implementation of photon mapping
class PhotonMapping : public Integrator {
 private:
  const int nPhotons;
  const int nDensityEstimation;
  const int maxDepth;
  PhotonMap photonMap;

 public:
  PhotonMapping(int nPhotons, int nDensityEstimation, int maxDepth = 100)
      : nPhotons(nPhotons),
        nDensityEstimation(nDensityEstimation),
        maxDepth(maxDepth) {}

  const PhotonMap* getPhotonMapPtr() const { return &photonMap; }

  void build(const Scene& scene, Sampler& sampler) override {
    std::vector<Photon> photons;

    // photon tracing
    spdlog::info("[PhotonMapping] tracing photons");
#pragma omp parallel for schedule(dynamic, 1)
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
      const Vec3 dir =
          light->sampleDirection(light_surf, sampler, light_dir_pdf);

      // spawn ray
      Ray ray(light_surf.position, dir);
      Vec3 throughput = light->Le(light_surf, dir) /
                        (light_choose_pdf * light_pos_pdf * light_dir_pdf) *
                        std::abs(dot(dir, light_surf.normal));

      // trace photons
      // whener hitting diffuse surface, add photon to the photon array
      // recursively tracing photon with russian roulette
      for (int k = 0; k < maxDepth; ++k) {
        IntersectInfo info;
        if (scene.intersect(ray, info)) {
          // if hitting diffuse surface, add photon to the photon array
          const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();
          if (bxdf_type == BxDFType::DIFFUSE) {
            // TODO: remove lock to get more speed
#pragma omp critical
            {
              photons.emplace_back(throughput, info.surfaceInfo.position,
                                   -ray.direction);
            }
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
          const Vec3 f = info.hitPrimitive->sampleBxDF(
              -ray.direction, info.surfaceInfo, sampler, dir, pdf_dir);

          // update throughput and ray
          throughput *=
              f * std::abs(dot(dir, info.surfaceInfo.normal)) / pdf_dir;
          ray = Ray(info.surfaceInfo.position, dir);
        } else {
          // photon goes to the sky
          break;
        }
      }
    }

    // build photon map
    spdlog::info("[PhotonMapping] building photon map");
    photonMap.setPhotons(photons);
    photonMap.build();
  }

  Vec3 integrate(const Ray& ray_in, const Scene& scene,
                 Sampler& sampler) const override {
    // recursively raytrace until hitting diffuse surface
    Ray ray = ray_in;
    Vec3 throughput(1, 1, 1);
    for (int k = 0; k < maxDepth; ++k) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // when directly hitting light, break
        if (info.hitPrimitive->hasAreaLight()) {
          return throughput *
                 info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
        }

        const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();

        // if hitting diffuse surface, query nearby photons and compute
        // reflected radiance
        if (bxdf_type == BxDFType::DIFFUSE) {
          // get nearby photons
          float max_dist2;
          const std::vector<int> photon_indices =
              photonMap.queryKNearestPhotons(info.surfaceInfo.position,
                                             nDensityEstimation, max_dist2);

          // compute reflected radiance with simple kernel
          Vec3 Lo;
          for (const int photon_idx : photon_indices) {
            const Photon& photon = photonMap.getIthPhoton(photon_idx);
            const Vec3 f = info.hitPrimitive->evaluateBxDF(
                -ray.direction, photon.wi, info.surfaceInfo);
            Lo += f * photon.throughput;
          }
          Lo /= (nPhotons * PI * max_dist2);

          return throughput * Lo;
        }
        // generate next ray
        else if (bxdf_type == BxDFType::SPECULAR) {
          // sample direction by BxDF
          Vec3 dir;
          float pdf_dir;
          const Vec3 f = info.hitPrimitive->sampleBxDF(
              -ray.direction, info.surfaceInfo, sampler, dir, pdf_dir);

          // update throughput and ray
          throughput *=
              f * std::abs(dot(dir, info.surfaceInfo.normal)) / pdf_dir;
          ray = Ray(info.surfaceInfo.position, dir);
        } else {
          spdlog::error("[PhotonMapping] invalid BxDF type");
          break;
        }
      } else {
        break;
      }
    }
    return Vec3(0);
  }
};

#endif