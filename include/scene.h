#ifndef _SCENE_H
#define _SCENE_H
#include <memory>
#include <vector>

#include "core.h"
#include "primitive.h"

// linear intersector
class Intersector {
 private:
  const Primitive* primitives;
  unsigned int nPrimitives;

 public:
  Intersector() {}

  void setPrimitives(const std::vector<Primitive>& primitives) {
    this->primitives = primitives.data();
    this->nPrimitives = primitives.size();
  }

  bool intersect(const Ray& ray, IntersectInfo& info) const {
    bool hit = false;
    for (unsigned int i = 0; i < nPrimitives; ++i) {
      const Primitive& p = primitives[i];
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

  void build() { intersector.setPrimitives(primitives); }

  bool intersect(const Ray& ray, IntersectInfo& info) {
    return intersector.intersect(ray, info);
  }
};

#endif