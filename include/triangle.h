#ifndef _SHAPE_H
#define _SHAPE_H
#include "core.h"
#include "sampler.h"

class Triangle {
 private:
  const float* vertices;
  const uint32_t* indices;
  const float* normals;
  const float* texcoords;

  const uint32_t faceID;
  Vec3f geometricNormal;
  float surfaceArea;

 public:
  Triangle(const float* vertices, const uint32_t* indices, const float* normals,
           const float* texcoords, uint32_t faceID)
      : vertices(vertices),
        indices(indices),
        normals(normals),
        texcoords(texcoords),
        faceID(faceID) {
    const Vec3ui vidx = getIndices();
    const Vec3f p1 = getVertexPosition(vidx[0]);
    const Vec3f p2 = getVertexPosition(vidx[1]);
    const Vec3f p3 = getVertexPosition(vidx[2]);

    // compute geometric normal
    geometricNormal = normalize(cross(p2 - p1, p3 - p1));

    // compute surface area
    surfaceArea = 0.5f * length(cross(p2 - p1, p3 - p1));
  }

  // return vertex position
  Vec3f getVertexPosition(uint32_t vertexID) const {
    return Vec3f(vertices[3 * vertexID + 0], vertices[3 * vertexID + 1],
                 vertices[3 * vertexID + 2]);
  }

  // return vertex normal
  Vec3f getVertexNormal(uint32_t vertexID) const {
    return Vec3f(normals[3 * vertexID + 0], normals[3 * vertexID + 1],
                 normals[3 * vertexID + 2]);
  }

  // return vertex texcoords
  Vec2f getVertexTexcoords(uint32_t vertexID) const {
    return Vec2f(texcoords[2 * vertexID + 0], texcoords[2 * vertexID + 1]);
  }

  // return vertex indices
  Vec3ui getIndices() const {
    return Vec3ui(indices[3 * faceID + 0], indices[3 * faceID + 1],
                  indices[3 * faceID + 2]);
  }

  // return geometric normal
  Vec3f getGeometricNormal() const { return geometricNormal; }

  // compute shading normal at given position
  Vec3f computeShadingNormal(const Vec2f& barycentric) const {
    const Vec3ui vidx = getIndices();
    const Vec3f n1 = getVertexNormal(vidx[0]);
    const Vec3f n2 = getVertexNormal(vidx[1]);
    const Vec3f n3 = getVertexNormal(vidx[2]);
    return n1 * (1.0f - barycentric[0] - barycentric[1]) + n2 * barycentric[0] +
           n3 * barycentric[1];
  }

  // compute texcoords at given position
  Vec2f getTexcoords(const Vec2f& barycentric) const {
    const Vec3ui vidx = getIndices();
    const Vec2f t1 = getVertexTexcoords(vidx[0]);
    const Vec2f t2 = getVertexTexcoords(vidx[1]);
    const Vec2f t3 = getVertexTexcoords(vidx[2]);
    return t1 * (1.0f - barycentric[0] - barycentric[1]) + t2 * barycentric[0] +
           t3 * barycentric[1];
  }

  // sample point on the triangle
  SurfaceInfo samplePoint(Sampler& sampler, float& pdf) const {
    SurfaceInfo ret;

    const Vec3ui vidx = getIndices();

    const Vec3f p1 = getVertexPosition(vidx[0]);
    const Vec3f p2 = getVertexPosition(vidx[1]);
    const Vec3f p3 = getVertexPosition(vidx[2]);

    // sample point on triangle
    const Vec2f uv = sampler.getNext2D();
    const float su0 = std::sqrt(uv[0]);
    const Vec2f barycentric = Vec2f(1.0f - su0, uv[1] * su0);
    ret.position = (1.0f - barycentric[0] - barycentric[1]) * p1 +
                   barycentric[0] * p2 + barycentric[1] * p3;

    // compute normal
    ret.shadingNormal = computeShadingNormal(barycentric);

    // compute dpdu, dpdv
    orthonormalBasis(ret.shadingNormal, ret.dpdu, ret.dpdv);

    // compute texcoords
    ret.texcoords = getTexcoords(barycentric);
    ret.barycentric = barycentric;

    // compute pdf
    pdf = 1.0f / surfaceArea;

    return ret;
  }
};

#endif