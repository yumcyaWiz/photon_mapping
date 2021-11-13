#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H
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
class PM : public Integrator {
 private:
  const int n_photons;
  const int n_density_estimation;
  static constexpr int max_depth = 100;
  PhotonMap photon_map;

 public:
  PM(int nPhotons, int nDensityEstimation)
      : n_photons(nPhotons), n_density_estimation(nDensityEstimation) {}

  void build(const Scene& scene, Sampler& sampler) override {
    // TODO: trace photons, build photon map
    std::vector<Photon> photons(n_photons);

    // photon tracing
    for (int i = 0; i < n_photons; ++i) {
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

      // spawn ray
      Ray ray(light_surf.position, light_dir_pdf);
      Vec3 throughput = light->Le(light_surf, light_dir) /
                        (light_choose_pdf * light_pos_pdf) *
                        dot(light_dir, light_surf.normal);

      // trace photons
      // whener hitting diffuse surface, add photon to the photon array
      // recursively tracing photon with russian roulette
      for (int k = 0; k < max_depth; ++k) {
        IntersectInfo info;
        if (scene.intersect(ray, info)) {
          // if hitting diffuse surface, add photon to the photon array
          BxDFType bxdf_type = info.hitPrimitive->getBxDFType();
          if (bxdf_type == BxDFType::DIFFUSE) {
            photons.emplace_back(throughput, info.surfaceInfo.position,
                                 -ray.direction);
          }

          // russian roulette
          const float throughput_max =
              std::max(throughput[0], std::max(throughput[1], throughput[2]));
          if (sampler.getNext1D() > throughput_max) {
            break;
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

    // build photon map
    photon_map.setPhotons(photons);
    photon_map.build();
  }

  Vec3 integrate(const Ray& ray, const Scene& scene,
                 Sampler& sampler) const override {
    // TODO: trace from camera until hitting diffuse surface, compute radiance
    return Vec3(0);
  }
};

#endif