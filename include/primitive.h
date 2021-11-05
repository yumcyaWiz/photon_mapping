#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H
#include <cmath>
#include <memory>

#include "core.h"
#include "material.h"

class Shape {
 public:
  virtual bool intersect(const Ray& ray, IntersectInfo& info) const = 0;
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
    info.hitPos = ray(t);
    info.hitNormal = normalize(info.hitPos - center);
    return true;
  }
};

class Plane : public Shape {
 private:
  const Vec3 leftCornerPoint;
  const Vec3 right;
  const Vec3 up;

 public:
  Plane(const Vec3& leftCornerPoint, const Vec3& right, const Vec3& up)
      : leftCornerPoint(leftCornerPoint), right(right), up(up) {}

  bool intersect(const Ray& ray, IntersectInfo& info) const override {
    const Vec3 normal = normalize(cross(right, up));
    const Vec3 center = leftCornerPoint + 0.5f * right + 0.5f * up;
    const Vec3 rightDir = normalize(right);
    const float rightLength = length(right);
    const Vec3 upDir = normalize(up);
    const float upLength = length(up);

    const float t =
        -dot(ray.origin - center, normal) / dot(ray.direction, normal);
    if (t < ray.tmin || t > ray.tmax) return false;

    const Vec3 hitPos = ray(t);
    const float dx = dot(hitPos - leftCornerPoint, rightDir);
    const float dy = dot(hitPos - leftCornerPoint, upDir);
    if (dx < 0.0f || dx > rightLength || dy < 0.0f || dy > upLength)
      return false;

    info.t = t;
    info.hitPos = hitPos;
    info.hitNormal = normal;
    return true;
  }
};

class Primitive {
 private:
  std::shared_ptr<Shape> shape;
  std::shared_ptr<Material> material;

 public:
  Primitive(const std::shared_ptr<Shape>& shape,
            const std::shared_ptr<Material>& material)
      : shape(shape), material(material) {}

  bool intersect(const Ray& ray, IntersectInfo& info) const {
    return shape->intersect(ray, info);
  }

  Vec3 sampleBRDF(const Vec3& wo, Vec3& wi, float& pdf);
};

#endif