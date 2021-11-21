#ifndef _SHAPE_H
#define _SHAPE_H
#include "core.h"
#include "sampler.h"

class Shape {
 public:
  virtual SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const = 0;
};

class Triangle : public Shape {
 private:
  const float* vertices;
  const uint32_t* indices;
  const float* normals;
  const float* texcoords;

  const uint32_t faceID;

 public:
  Triangle(const float* vertices, const uint32_t* indices, const float* normals,
           const float* texcoords, uint32_t faceID)
      : vertices(vertices),
        indices(indices),
        normals(normals),
        texcoords(texcoords),
        faceID(faceID) {}

  SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const override {}
};

#endif