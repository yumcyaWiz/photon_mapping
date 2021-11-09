#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H
#include <cmath>
#include <memory>

#include "core.h"
#include "material.h"
#include "sampler.h"

class Shape {
 public:
  virtual bool intersect(const Ray& ray, IntersectInfo& info) const = 0;
  virtual SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const = 0;
};

class Sphere : public Shape {
 private:
  const Vec3 center;
  const float radius;

 public:
  Sphere(const Vec3& center, float radius) : center(center), radius(radius) {}

  bool intersect(const Ray& ray, IntersectInfo& info) const override {
    const float b = dot(ray.origin - center, ray.direction);
    const float c = length2(ray.origin - center) - radius * radius;
    const float D = b * b - c;
    if (D < 0) return false;

    const float t1 = -b - std::sqrt(D);
    const float t2 = -b + std::sqrt(D);
    float t = t1;
    if (t < ray.tmin || t > ray.tmax) {
      t = t2;
      if (t < ray.tmin || t > ray.tmax) return false;
    }

    info.t = t;
    info.surfaceInfo.position = ray(t);
    info.surfaceInfo.normal = normalize(info.surfaceInfo.position - center);
    return true;
  }

  SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const override {
    SurfaceInfo ret;
    const Vec3 p = sampleSphere(sampler.getNext2D(), pdf);
    ret.position = center + radius * p;
    ret.normal = p;
    pdf /= (radius * radius);
    return ret;
  }
};

class Plane : public Shape {
 private:
  const Vec3 leftCornerPoint;
  Vec3 rightDir;
  float rightLength;
  Vec3 upDir;
  float upLength;
  Vec3 center;
  Vec3 normal;

 public:
  Plane(const Vec3& leftCornerPoint, const Vec3& right, const Vec3& up)
      : leftCornerPoint(leftCornerPoint) {
    rightDir = normalize(right);
    rightLength = length(right);
    upDir = normalize(up);
    upLength = length(up);
    center = leftCornerPoint + 0.5f * rightLength * rightDir +
             0.5f * upLength * upDir;
    normal = normalize(cross(rightDir, upDir));
  }

  bool intersect(const Ray& ray, IntersectInfo& info) const override {
    const float t =
        -dot(ray.origin - center, normal) / dot(ray.direction, normal);
    if (t < ray.tmin || t > ray.tmax) return false;

    const Vec3 hitPos = ray(t);
    const float dx = dot(hitPos - leftCornerPoint, rightDir);
    const float dy = dot(hitPos - leftCornerPoint, upDir);
    if (dx < 0.0f || dx > rightLength || dy < 0.0f || dy > upLength)
      return false;

    info.t = t;
    info.surfaceInfo.position = hitPos;
    info.surfaceInfo.normal = normal;
    return true;
  }

  SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const override {
    SurfaceInfo ret;
    const Vec2 p = samplePlane(sampler.getNext2D(), rightLength, upLength, pdf);
    ret.position = leftCornerPoint + p[0] * rightDir + p[1] * upDir;
    ret.normal = normal;
    return ret;
  }
};

class Primitive {
 private:
  const std::shared_ptr<Shape> shape;
  const std::shared_ptr<Material> material;

 public:
  Primitive(const std::shared_ptr<Shape>& shape,
            const std::shared_ptr<Material>& material)
      : shape(shape), material(material) {}

  bool intersect(const Ray& ray, IntersectInfo& info) const {
    if (shape->intersect(ray, info)) {
      info.hitPrimitive = this;
      return true;
    }
    return false;
  }

  Vec3 sampleBRDF(const Vec3& wo, Vec3& wi, float& pdf) const;
};

#endif