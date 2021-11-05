#ifndef _SCENE_H
#define _SCENE_H
#include <memory>
#include <vector>

#include "core.h"
#include "primitive.h"

// linear intersector
class Intersector {
 private:
  std::shared_ptr<std::vector<Primitive>> primitives;

 public:
  Intersector() {}

  void setPrimitives(
      const std::shared_ptr<std::vector<Primitive>>& primitives) {
    this->primitives = primitives;
  }

  bool intersect(const Ray& ray, IntersectInfo& info) const {
    bool hit = false;
    for (const auto& p : *primitives) {
      if (p.intersect(ray, info)) {
        hit = true;
        ray.tmax = info.t;
      }
    }
    return hit;
  }
};

class Scene {
 private:
  std::vector<Primitive> primitives;
  Intersector intersector;

 public:
  Scene() {}

  void addPrimitive(const Primitive& primitive) {
    primitives.push_back(primitive);
  }

  void build() {
    intersector.setPrimitives(
        std::make_shared<std::vector<Primitive>>(this->primitives));
  }

  bool intersect(const Ray& ray, IntersectInfo& info) {
    return intersector.intersect(ray, info);
  }
};

#endif