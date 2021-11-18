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
// NOTE: for reference purpose
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
  const int nPhotonsGlobal;
  const int nEstimationGlobal;
  const int nPhotonsCaustics;
  const int nEstimationCaustics;
  const int maxDepth;
  bool finalGathering;

  PhotonMap globalPhotonMap;
  PhotonMap causticsPhotonMap;

  // compute reflected radiance with global photon map
  Vec3 computeRadianceWithPhotonMap(const Vec3& wo,
                                    const IntersectInfo& info) const {
    // get nearby photons
    float max_dist2;
    const std::vector<int> photon_indices =
        globalPhotonMap.queryKNearestPhotons(info.surfaceInfo.position,
                                             nEstimationGlobal, max_dist2);

    Vec3 Lo;
    for (const int photon_idx : photon_indices) {
      const Photon& photon = globalPhotonMap.getIthPhoton(photon_idx);
      const Vec3 f =
          info.hitPrimitive->evaluateBxDF(wo, photon.wi, info.surfaceInfo);
      Lo += f * photon.throughput;
    }
    if (photon_indices.size() > 0) {
      Lo /= (nPhotonsGlobal * PI * max_dist2);
    }
    return Lo;
  }

  // compute reflected radiance with caustics photon map
  Vec3 computeCausticsWithPhotonMap(const Vec3& wo,
                                    const IntersectInfo& info) const {
    // get nearby photons
    float max_dist2;
    const std::vector<int> photon_indices =
        causticsPhotonMap.queryKNearestPhotons(info.surfaceInfo.position,
                                               nEstimationGlobal, max_dist2);

    Vec3 Lo;
    for (const int photon_idx : photon_indices) {
      const Photon& photon = causticsPhotonMap.getIthPhoton(photon_idx);
      const Vec3 f =
          info.hitPrimitive->evaluateBxDF(wo, photon.wi, info.surfaceInfo);
      Lo += f * photon.throughput;
    }
    if (photon_indices.size() > 0) {
      Lo /= (nPhotonsCaustics * PI * max_dist2);
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
    Ray ray_shadow(info.surfaceInfo.position, wi);
    ray_shadow.tmax = r - RAY_EPS;

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

  Vec3 computeIndirectIlluminationRecursive(const Scene& scene, const Vec3& wo,
                                            const IntersectInfo& info,
                                            Sampler& sampler, int depth) const {
    if (depth >= maxDepth) return Vec3(0);

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
      const BxDFType bxdf_type = info_fg.hitPrimitive->getBxDFType();

      // when hitting diffuse, compute radiance with photon map
      if (bxdf_type == BxDFType::DIFFUSE) {
        Li += f * cos *
              computeRadianceWithPhotonMap(-ray_fg.direction, info_fg) /
              pdf_dir;
      }
      // when hitting specular, recursively call this function
      // NOTE: to include the path like LSDSDE
      else if (bxdf_type == BxDFType::SPECULAR) {
        Li += f * cos *
              computeIndirectIlluminationRecursive(
                  scene, -ray_fg.direction, info_fg, sampler, depth + 1) /
              pdf_dir;
      }
    }

    return Li;
  }

  // compute indirect illumination with final gathering
  Vec3 computeIndirectIllumination(const Scene& scene, const Vec3& wo,
                                   const IntersectInfo& info,
                                   Sampler& sampler) const {
    return computeIndirectIlluminationRecursive(scene, wo, info, sampler, 0);
  }

  // sample initial ray from light and compute initial throughput
  Ray sampleRayFromLight(const Scene& scene, Sampler& sampler,
                         Vec3& throughput) {
    // sample light
    float light_choose_pdf;
    const std::shared_ptr<Light> light =
        scene.sampleLight(sampler, light_choose_pdf);

    // sample point on light
    float light_pos_pdf;
    const SurfaceInfo light_surf = light->samplePoint(sampler, light_pos_pdf);

    // sample direction on light
    float light_dir_pdf;
    const Vec3 dir = light->sampleDirection(light_surf, sampler, light_dir_pdf);

    // spawn ray
    Ray ray(light_surf.position, dir);
    throughput = light->Le(light_surf, dir) /
                 (light_choose_pdf * light_pos_pdf * light_dir_pdf) *
                 std::abs(dot(dir, light_surf.normal));

    return ray;
  }

 public:
  PhotonMapping(int nPhotonsGlobal, int nEstimationGlobal,
                float nPhotonsCausticsMultiplier, int nEstimationCaustics,
                bool finalGathering, int maxDepth)
      : nPhotonsGlobal(nPhotonsGlobal),
        nEstimationGlobal(nEstimationGlobal),
        nPhotonsCaustics(nPhotonsGlobal * nPhotonsCausticsMultiplier),
        nEstimationCaustics(nEstimationCaustics),
        finalGathering(finalGathering),
        maxDepth(maxDepth) {}

  const PhotonMap* getPhotonMapPtr() const { return &globalPhotonMap; }

  // photon tracing and build photon map
  void build(const Scene& scene, Sampler& sampler) override {
    std::vector<Photon> photons;

    // init sampler for each thread
    std::vector<std::unique_ptr<Sampler>> samplers(omp_get_max_threads());
    for (int i = 0; i < samplers.size(); ++i) {
      samplers[i] = sampler.clone();
      samplers[i]->setSeed(samplers[i]->getSeed() * (i + 1));
    }

    // build global photon map
    // photon tracing
    spdlog::info(
        "[PhotonMapping] tracing photons for building global photon map");
#pragma omp parallel for
    for (int i = 0; i < nPhotonsGlobal; ++i) {
      auto& sampler_per_thread = *samplers[omp_get_thread_num()];

      // sample initial ray from light and set initial throughput
      Vec3 throughput;
      Ray ray = sampleRayFromLight(scene, sampler_per_thread, throughput);

      // trace photons
      // whener hitting diffuse surface, add photon to the photon array
      // recursively tracing photon with russian roulette
      for (int k = 0; k < maxDepth; ++k) {
        IntersectInfo info;
        if (scene.intersect(ray, info)) {
          const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();
          if (bxdf_type == BxDFType::DIFFUSE &&
              !info.hitPrimitive->hasAreaLight()) {
            if (std::isnan(throughput[0]) || std::isnan(throughput[1]) ||
                std::isnan(throughput[2])) {
              spdlog::error("[PhotonMapping] photon throughput is NaN");
              break;
            }

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
    spdlog::info("[PhotonMapping] building global photon map");
    globalPhotonMap.setPhotons(photons);
    globalPhotonMap.build();

    // build caustics photon map
    if (finalGathering) {
      photons.clear();

      // photon tracing
      spdlog::info(
          "[PhotonMapping] tracing photons for building caustics photon map");
#pragma omp parallel for
      for (int i = 0; i < nPhotonsCaustics; ++i) {
        auto& sampler_per_thread = *samplers[omp_get_thread_num()];

        // sample initial ray from light and set initial throughput
        Vec3 throughput;
        Ray ray = sampleRayFromLight(scene, sampler_per_thread, throughput);

        // when hitting diffuse surface after specular, add photon to the photon
        // array
        bool prev_specular = false;
        for (int k = 0; k < maxDepth; ++k) {
          IntersectInfo info;
          if (scene.intersect(ray, info)) {
            const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();
            if (prev_specular && bxdf_type == BxDFType::DIFFUSE &&
                !info.hitPrimitive->hasAreaLight()) {
              if (std::isnan(throughput[0]) || std::isnan(throughput[1]) ||
                  std::isnan(throughput[2])) {
                spdlog::error("[PhotonMapping] photon throughput is NaN");
                break;
              }

              // TODO: remove lock to get more speed
#pragma omp critical
              {
                photons.emplace_back(throughput, info.surfaceInfo.position,
                                     -ray.direction);
              }

              // break after the first diffuse surface
              break;
            }

            prev_specular = (bxdf_type == BxDFType::SPECULAR);

            // russian roulette
            if (k > 0) {
              const float russian_roulette_prob =
                  std::min(std::max(throughput[0],
                                    std::max(throughput[1], throughput[2])),
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

      spdlog::info("[PhotonMapping] building caustics photon map");
      causticsPhotonMap.setPhotons(photons);
      causticsPhotonMap.build();
    }
  }

  Vec3 integrate(const Ray& ray_in, const Scene& scene,
                 Sampler& sampler) const override {
    Ray ray = ray_in;
    Vec3 throughput(1, 1, 1);

    // recursively raytrace until hitting diffuse surface
    for (int k = 0; k < maxDepth; ++k) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // when directly hitting light, break
        if (info.hitPrimitive->hasAreaLight()) {
          return throughput *
                 info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
        }

        const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();

        // if hitting diffuse surface, computed reflected radiance with photon
        // map
        if (bxdf_type == BxDFType::DIFFUSE) {
          if (!finalGathering) {
            return throughput *
                   computeRadianceWithPhotonMap(-ray.direction, info);
          } else {
            // compute direct illumination by explicit light sampling
            const Vec3 Ld =
                computeDirectIllumination(scene, -ray.direction, info, sampler);

            // compute caustics illumination with caustics photon map
            const Vec3 Lc = computeCausticsWithPhotonMap(-ray.direction, info);

            // compute indirect illumination with final gathering
            const Vec3 Li = computeIndirectIllumination(scene, -ray.direction,
                                                        info, sampler);

            return throughput * (Ld + Lc + Li);
          }
        }
        // if hitting specular surface, generate next ray and continue
        // raytracing
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