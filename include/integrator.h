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
  bool finalGathering;

  PhotonMap photonMap;

  // compute reflected radiance with photon map
  Vec3 computeReflectedRadiance(const Vec3& wo, const IntersectInfo& info,
                                const std::vector<int>& photon_indices,
                                float max_dist2) const {
    Vec3 Lo;
    for (const int photon_idx : photon_indices) {
      const Photon& photon = photonMap.getIthPhoton(photon_idx);
      const Vec3 f =
          info.hitPrimitive->evaluateBxDF(wo, photon.wi, info.surfaceInfo);
      Lo += f * photon.throughput;
    }
    if (photon_indices.size() > 0) {
      Lo /= (nPhotons * PI * max_dist2);
    }
    return Lo;
  }

  // compute direct illumination with explicit light sampling(NEE)
  Vec3 computeDirectIllumination(const Scene& scene, const Vec3& wo,
                                 const IntersectInfo& info,
                                 Sampler& sampler) const {
    Vec3 Ld;

    // sample light
    float pdf_choose_light;
    const std::shared_ptr<Light> light =
        scene.sampleLight(sampler, pdf_choose_light);

    // sample point on light
    float pdf_pos_light;
    const SurfaceInfo light_surf = light->samplePoint(sampler, pdf_pos_light);

    // convert positional pdf to directional pdf
    const Vec3 wi = normalize(light_surf.position - info.surfaceInfo.position);
    const float r = length(light_surf.position - info.surfaceInfo.position);
    const float pdf_dir =
        pdf_pos_light * r * r / std::abs(dot(-wi, light_surf.normal));

    // create shadow ray
    constexpr float EPS = 0.001f;
    Ray ray_shadow(info.surfaceInfo.position, wi);
    ray_shadow.tmax = r - EPS;

    // trace ray to the light
    IntersectInfo info_shadow;
    if (!scene.intersect(ray_shadow, info_shadow)) {
      const Vec3 Le = light->Le(light_surf, -wi);
      const Vec3 f = info.hitPrimitive->evaluateBxDF(wo, wi, info.surfaceInfo);
      const float cos = std::abs(dot(wi, info.surfaceInfo.normal));
      Ld = f * cos * Le / (pdf_choose_light * pdf_dir);
    }

    return Ld;
  }

  // compute indirect illumination with final gathering
  Vec3 computeIndirectIllumination(const Scene& scene, const Vec3& wo,
                                   const IntersectInfo& info,
                                   Sampler& sampler) const {
    Vec3 Li;

    // sample direction by BxDF
    Vec3 dir;
    float pdf_dir;
    const Vec3 f = info.hitPrimitive->sampleBxDF(wo, info.surfaceInfo, sampler,
                                                 dir, pdf_dir);
    const float cos = std::abs(dot(info.surfaceInfo.normal, dir));

    // trace final gathering ray
    Ray ray_fg(info.surfaceInfo.position, dir);
    IntersectInfo info_fg;
    if (scene.intersect(ray_fg, info_fg)) {
      if (info_fg.hitPrimitive->getBxDFType() == BxDFType::DIFFUSE) {
        // get nearby photons
        float max_dist2;
        const std::vector<int> photon_indices = photonMap.queryKNearestPhotons(
            info_fg.surfaceInfo.position, nDensityEstimation, max_dist2);

        Li += f * cos *
              computeReflectedRadiance(-ray_fg.direction, info_fg,
                                       photon_indices, max_dist2) /
              pdf_dir;
      }
    }

    return Li;
  }

 public:
  PhotonMapping(int nPhotons, int nDensityEstimation, bool finalGathering,
                int maxDepth = 100)
      : nPhotons(nPhotons),
        nDensityEstimation(nDensityEstimation),
        finalGathering(finalGathering),
        maxDepth(maxDepth) {}

  const PhotonMap* getPhotonMapPtr() const { return &photonMap; }

  void build(const Scene& scene, Sampler& sampler) override {
    std::vector<Photon> photons;

    // init sampler for each thread
    std::vector<std::unique_ptr<Sampler>> samplers(omp_get_max_threads());
    for (int i = 0; i < samplers.size(); ++i) {
      samplers[i] = sampler.clone();
      samplers[i]->setSeed(samplers[i]->getSeed() * (i + 1));
    }

    // photon tracing
    spdlog::info("[PhotonMapping] tracing photons");
#pragma omp parallel for
    for (int i = 0; i < nPhotons; ++i) {
      auto& sampler_per_thread = *samplers[omp_get_thread_num()];

      // sample light
      float light_choose_pdf;
      const std::shared_ptr<Light> light =
          scene.sampleLight(sampler_per_thread, light_choose_pdf);

      // sample point on light
      float light_pos_pdf;
      const SurfaceInfo light_surf =
          light->samplePoint(sampler_per_thread, light_pos_pdf);

      // sample direction on light
      float light_dir_pdf;
      const Vec3 dir =
          light->sampleDirection(light_surf, sampler_per_thread, light_dir_pdf);

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
            if (sampler_per_thread.getNext1D() >= russian_roulette_prob) {
              break;
            }
            throughput /= russian_roulette_prob;
          }

          // sample direction by BxDF
          Vec3 dir;
          float pdf_dir;
          const Vec3 f =
              info.hitPrimitive->sampleBxDF(-ray.direction, info.surfaceInfo,
                                            sampler_per_thread, dir, pdf_dir);

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
          if (!finalGathering) {
            // get nearby photons
            float max_dist2;
            const std::vector<int> photon_indices =
                photonMap.queryKNearestPhotons(info.surfaceInfo.position,
                                               nDensityEstimation, max_dist2);

            return throughput * computeReflectedRadiance(-ray.direction, info,
                                                         photon_indices,
                                                         max_dist2);
          }
          // final gathering
          // trace one more ray, and compute reflected radiance there
          else {
            // compute direct illumination
            Vec3 Ld =
                computeDirectIllumination(scene, -ray.direction, info, sampler);

            // compute indirect illumination with final gathering
            Vec3 Li = computeIndirectIllumination(scene, -ray.direction, info,
                                                  sampler);

            return throughput * (Ld + Li);
          }
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