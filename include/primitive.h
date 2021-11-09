#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H
#include <cmath>
#include <memory>

#include "material.h"
#include "shape.h"

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