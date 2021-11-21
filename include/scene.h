#ifndef _SCENE_H
#define _SCENE_H
#include <embree3/rtcore.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <memory>
#include <vector>

#include "core.h"
#include "primitive.h"
#include "tiny_obj_loader.h"

// create material from tinyobj material
Material createMaterial(const tinyobj::material_t& material) {
  const Vec3f kd =
      Vec3f(material.diffuse[0], material.diffuse[1], material.diffuse[2]);
  const Vec3f ks =
      Vec3f(material.specular[0], material.specular[1], material.specular[2]);
  const Vec3f ke =
      Vec3f(material.emission[0], material.emission[1], material.emission[2]);
  return Material(kd, ks, ke);
}

class Scene {
 private:
  // triangles
  // NOTE: size of normals, texcoords == size of vertices
  std::vector<float> vertices;
  std::vector<uint32_t> indices;
  std::vector<float> normals;
  std::vector<float> texcoords;

  // materials
  std::vector<Material> materials;

  // lights
  std::vector<std::shared_ptr<Light>> lights;

  // embree
  RTCDevice device;
  RTCScene scene;

  // return vertex normal
  Vec3f getVertexNormal(uint32_t vertexID) const {
    return Vec3f(normals[3 * vertexID + 0], normals[3 * vertexID + 1],
                 normals[3 * vertexID + 2]);
  }

  // return vertex texcoords
  Vec2f getVertexTexcoords(uint32_t vertexID) const {
    return Vec2f(texcoords[2 * vertexID + 0], texcoords[2 * vertexID + 1]);
  }

  // return vertex indices of specified face
  Vec3i getIndices(uint32_t faceID) const {
    return Vec3i(indices[3 * faceID + 0], indices[3 * faceID + 1],
                 indices[3 * faceID + 2]);
  }

  // compute normal of specified face, barycentric
  Vec3f getFaceNormal(uint32_t faceID, const Vec2f& barycentric) const {
    const Vec3i vidx = getIndices(faceID);
    const Vec3f n1 = getVertexNormal(vidx[0]);
    const Vec3f n2 = getVertexNormal(vidx[1]);
    const Vec3f n3 = getVertexNormal(vidx[2]);
    return n1 * (1.0f - barycentric[0] - barycentric[1]) + n2 * barycentric[0] +
           n3 * barycentric[1];
  }

  // compute texcoords of specified face, barycentric
  Vec2f getTexcoords(uint32_t faceID, const Vec2f& barycentric) const {
    const Vec3i vidx = getIndices(faceID);
    const Vec2f t1 = getVertexTexcoords(vidx[0]);
    const Vec2f t2 = getVertexTexcoords(vidx[1]);
    const Vec2f t3 = getVertexTexcoords(vidx[2]);
    return t1 * (1.0f - barycentric[0] - barycentric[1]) + t2 * barycentric[0] +
           t3 * barycentric[1];
  }

  // get material of specified face
  const Material* getMaterial(uint32_t faceID) const {
    if (materials.size() > 0) {
      return &materials[faceID];
    }
  }

  void clear() {
    vertices.clear();
    indices.clear();
    normals.clear();
    texcoords.clear();
    materials.clear();
  }

 public:
  Scene() {}

  // load obj file
  // TODO: remove vertex duplication
  void loadModel(const std::filesystem::path& filepath) {
    clear();

    spdlog::info("[Scene] loading: {}", filepath.generic_string());

    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = "./";
    reader_config.triangulate = true;

    tinyobj::ObjReader reader;
    if (!reader.ParseFromFile(filepath.generic_string(), reader_config)) {
      if (!reader.Error().empty()) {
        spdlog::error("[Scene] failed to load {} : {}",
                      filepath.generic_string(), reader.Error());
      }
      return;
    }

    if (!reader.Warning().empty()) {
      spdlog::warn("[Scene] {}", reader.Warning());
    }

    const auto& attrib = reader.GetAttrib();
    const auto& shapes = reader.GetShapes();
    const auto& materials = reader.GetMaterials();

    // loop over shapes
    for (size_t s = 0; s < shapes.size(); ++s) {
      size_t index_offset = 0;
      // loop over faces
      for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); ++f) {
        const size_t fv =
            static_cast<size_t>(shapes[s].mesh.num_face_vertices[f]);

        std::vector<Vec3f> vertices;
        std::vector<Vec3f> normals;
        std::vector<Vec2f> texcoords;

        // loop over vertices
        // get vertices, normals, texcoords of a triangle
        for (size_t v = 0; v < fv; ++v) {
          const tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

          const tinyobj::real_t vx =
              attrib.vertices[3 * static_cast<size_t>(idx.vertex_index) + 0];
          const tinyobj::real_t vy =
              attrib.vertices[3 * static_cast<size_t>(idx.vertex_index) + 1];
          const tinyobj::real_t vz =
              attrib.vertices[3 * static_cast<size_t>(idx.vertex_index) + 2];
          vertices.push_back(Vec3f(vx, vy, vz));

          if (idx.normal_index >= 0) {
            const tinyobj::real_t nx =
                attrib.normals[3 * static_cast<size_t>(idx.normal_index) + 0];
            const tinyobj::real_t ny =
                attrib.normals[3 * static_cast<size_t>(idx.normal_index) + 1];
            const tinyobj::real_t nz =
                attrib.normals[3 * static_cast<size_t>(idx.normal_index) + 2];
            normals.push_back(Vec3f(nx, ny, nz));
          }

          if (idx.texcoord_index >= 0) {
            const tinyobj::real_t tx =
                attrib
                    .texcoords[2 * static_cast<size_t>(idx.texcoord_index) + 0];
            const tinyobj::real_t ty =
                attrib
                    .texcoords[2 * static_cast<size_t>(idx.texcoord_index) + 1];
            texcoords.push_back(Vec2f(tx, ty));
          }
        }

        // if normals is empty, add geometric normal
        if (normals.size() == 0) {
          const Vec3f v1 = normalize(vertices[1] - vertices[0]);
          const Vec3f v2 = normalize(vertices[2] - vertices[0]);
          const Vec3f n = normalize(cross(v1, v2));
          normals.push_back(n);
          normals.push_back(n);
          normals.push_back(n);
        }

        // if texcoords is empty, add barycentric coords
        if (texcoords.size() == 0) {
          texcoords.push_back(Vec2f(0, 0));
          texcoords.push_back(Vec2f(1, 0));
          texcoords.push_back(Vec2f(0, 1));
        }

        for (int i = 0; i < 3; ++i) {
          this->vertices.push_back(vertices[i][0]);
          this->vertices.push_back(vertices[i][1]);
          this->vertices.push_back(vertices[i][2]);

          this->normals.push_back(normals[i][0]);
          this->normals.push_back(normals[i][1]);
          this->normals.push_back(normals[i][2]);

          this->texcoords.push_back(texcoords[i][0]);
          this->texcoords.push_back(texcoords[i][1]);

          this->indices.push_back(this->indices.size());
        }

        // add material
        // TODO: remove duplicate
        const int materialID = shapes[s].mesh.material_ids[f];
        if (materialID != -1) {
          const tinyobj::material_t& m = materials[materialID];
          this->materials.push_back(createMaterial(m));
        }

        index_offset += fv;
      }
    }

    spdlog::info("[Scene] vertices: {}", nVertices());
    spdlog::info("[Scene] faces: {}", nFaces());
    spdlog::info("[Scene] materials: {}", nMaterials());
  }

  uint32_t nVertices() const { return vertices.size() / 3; }
  uint32_t nFaces() const { return indices.size() / 3; }
  uint32_t nMaterials() const { return materials.size(); }

  void build() {
    spdlog::info("[Scene] building scene...", lights.size());

    // setup embree
    device = rtcNewDevice(NULL);
    scene = rtcNewScene(device);

    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    // set vertices
    float* vb = (float*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0,
                                                RTC_FORMAT_FLOAT3,
                                                3 * sizeof(float), nVertices());
    for (int i = 0; i < vertices.size(); ++i) {
      vb[i] = vertices[i];
    }

    // set indices
    unsigned* ib = (unsigned*)rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned),
        nFaces());
    for (int i = 0; i < indices.size(); ++i) {
      ib[i] = indices[i];
    }

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);

    // populate lights
    for (const auto& primitive : primitives) {
      if (primitive.hasAreaLight()) {
        lights.push_back(primitive.getAreaLightPtr());
      }
    }

    spdlog::info("[Scene] number of lights: {}", lights.size());
  }

  bool intersect(const Ray& ray, IntersectInfo& info) const {
    RTCRayHit rayhit;
    rayhit.ray.org_x = ray.origin[0];
    rayhit.ray.org_y = ray.origin[1];
    rayhit.ray.org_z = ray.origin[2];
    rayhit.ray.dir_x = ray.direction[0];
    rayhit.ray.dir_y = ray.direction[1];
    rayhit.ray.dir_z = ray.direction[2];
    rayhit.ray.tnear = ray.tmin;
    rayhit.ray.tfar = ray.tmax;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    rtcIntersect1(scene, &context, &rayhit);

    if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
      info.t = rayhit.ray.tfar;
      info.primID = rayhit.hit.primID;

      // set surface info
      info.surfaceInfo.position = ray(info.t);
      info.surfaceInfo.barycentric = Vec2f(rayhit.hit.u, rayhit.hit.v);
      info.surfaceInfo.texcoords =
          getTexcoords(rayhit.hit.primID, info.surfaceInfo.barycentric);
      info.surfaceInfo.normal =
          getFaceNormal(rayhit.hit.primID, info.surfaceInfo.barycentric);
      orthonormalBasis(info.surfaceInfo.normal, info.surfaceInfo.dpdu,
                       info.surfaceInfo.dpdv);

      // set material
      info.material = getMaterial(rayhit.hit.primID);

      return true;
    } else {
      return false;
    }
  }

  std::shared_ptr<Light> sampleLight(Sampler& sampler, float& pdf) const {
    unsigned int lightIdx = lights.size() * sampler.getNext1D();
    if (lightIdx == lights.size()) lightIdx--;
    pdf = 1.0f / lights.size();
    return lights[lightIdx];
  }
};

#endif