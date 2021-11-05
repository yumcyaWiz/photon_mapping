#ifndef _CAMERA_H
#define _CAMERA_H
#include <cmath>

#include "spdlog/spdlog.h"
//
#include "core.h"

// pinhole camera
class Camera {
 private:
  Vec3 position;
  Vec3 forward;
  Vec3 right;
  Vec3 up;

  float FOV;
  float focal_length;

 public:
  Camera(const Vec3& position, const Vec3& forward, float FOV = 90.0f)
      : position(position), forward(forward) {
    right = normalize(cross(forward, Vec3(0, 1, 0)));
    up = normalize(cross(right, forward));

    spdlog::info("[Camera] position: ({}, {}, {})", position.x, position.y,
                 position.z);
    spdlog::info("[Camera] forward: ({}, {}, {})", forward.x, forward.y,
                 forward.z);
    spdlog::info("[Camera] right: ({}, {}, {})", right.x, right.y, right.z);
    spdlog::info("[Camera] up: ({}, {}, {})", up.x, up.y, up.z);

    // compute focal length from FOV
    focal_length = 1.0f / std::tan(deg2rad(0.5f * FOV));

    spdlog::info("[Camera] focal_length: {}", focal_length);
  }

  bool sampleRay(float u, float v, Ray& ray, float& pdf) const {
    const Vec3 pinholePos = position + focal_length * forward;
    const Vec3 sensorPos = position + u * right + v * up;
    ray = Ray(sensorPos, normalize(pinholePos - sensorPos));
    pdf = 1.0f;
    return true;
  }
};

#endif