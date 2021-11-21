#ifndef _SAMPLER_H
#define _SAMPLER_H
#include <cstdint>
#include <limits>
#include <memory>

#include "core.h"

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
typedef struct {
  uint64_t state;
  uint64_t inc;
} pcg32_random_t;

inline uint32_t pcg32_random_r(pcg32_random_t* rng) {
  uint64_t oldstate = rng->state;
  // Advance internal state
  rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
  // Calculate output function (XSH RR), uses old state for max ILP
  uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  uint32_t rot = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// random number generator
class RNG {
 private:
  pcg32_random_t state;

 public:
  RNG() {
    state.state = 1;
    state.inc = 1;
  }
  RNG(uint64_t seed) {
    state.state = seed;
    state.inc = 1;
  }

  uint64_t getSeed() const { return state.state; }
  void setSeed(uint64_t seed) { state.state = seed; }

  float getNext() {
    constexpr float divider = 1.0f / std::numeric_limits<uint32_t>::max();
    return pcg32_random_r(&state) * divider;
  }
};

// sampler interface
class Sampler {
 protected:
  RNG rng;

 public:
  Sampler() {}

  Sampler(uint64_t seed) : rng(seed) {}

  uint64_t getSeed() const { return rng.getSeed(); }
  void setSeed(uint64_t seed) { rng.setSeed(seed); }

  virtual std::unique_ptr<Sampler> clone() const = 0;
  virtual float getNext1D() = 0;
  virtual Vec2f getNext2D() = 0;
};

// uniform distribution sampler
class UniformSampler : public Sampler {
 public:
  UniformSampler() : Sampler() {}
  UniformSampler(uint64_t seed) : Sampler(seed) {}

  std::unique_ptr<Sampler> clone() const override {
    return std::make_unique<UniformSampler>();
  }

  float getNext1D() override { return rng.getNext(); }
  Vec2f getNext2D() override { return Vec2f(rng.getNext(), rng.getNext()); }
};

inline Vec3 sampleCosineHemisphere(const Vec2f& uv, float& pdf) {
  const float theta =
      0.5f * std::acos(std::clamp(1.0f - 2.0f * uv[0], -1.0f, 1.0f));
  const float phi = PI_MUL_2 * uv[1];
  const float cosTheta = std::cos(theta);
  pdf = PI_INV * cosTheta;
  return sphericalToCartesian(theta, phi);
}

inline Vec3 sampleSphere(const Vec2f& uv, float& pdf) {
  const float theta = std::acos(std::clamp(1.0f - 2.0f * uv[0], -1.0f, 1.0f));
  const float phi = PI_MUL_2 * uv[1];
  pdf = PI_MUL_4_INV;
  return sphericalToCartesian(theta, phi);
}

inline Vec2f samplePlane(const Vec2f& uv, float lx, float ly, float& pdf) {
  pdf = 1.0f / (lx * ly);
  return Vec2f(uv[0] * lx, uv[1] * ly);
}

#endif