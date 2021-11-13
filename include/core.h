#ifndef _CORE_H
#define _CORE_H
#include <cmath>
#include <iostream>
#include <limits>

constexpr float PI = 3.14159265359;

constexpr float PI_MUL_2 = 2.0f * PI;
constexpr float PI_MUL_4 = 4.0f * PI;

constexpr float PI_DIV_2 = 0.5f * PI;
constexpr float PI_DIV_4 = 0.25f * PI;

constexpr float PI_INV = 1.0f / PI;
constexpr float PI_MUL_2_INV = 1.0f / PI_MUL_2;
constexpr float PI_MUL_4_INV = 1.0f / PI_MUL_4;

inline float rad2deg(float rad) { return 180.0f * rad / PI; }
inline float deg2rad(float deg) { return deg / 180.0f * PI; }

struct Vec2 {
  float v[2];

  Vec2() { v[0] = v[1] = 0; }
  Vec2(float x) { v[0] = v[1] = x; }
  Vec2(float x, float y) {
    v[0] = x;
    v[1] = y;
  }

  float operator[](int i) const { return v[i]; }
  float& operator[](int i) { return v[i]; }
};

struct Vec3 {
  float v[3];

  Vec3() { v[0] = v[1] = v[2] = 0; }
  Vec3(float x) { v[0] = v[1] = v[2] = x; }
  Vec3(float x, float y, float z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }

  float operator[](int i) const { return v[i]; }
  float& operator[](int i) { return v[i]; }

  Vec3 operator-() const { return Vec3(-v[0], -v[1], -v[2]); }

  Vec3& operator*=(const Vec3& v) {
    this->v[0] *= v[0];
    this->v[1] *= v[1];
    this->v[2] *= v[2];
    return *this;
  }
};

inline Vec3 operator+(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
}
inline Vec3 operator+(const Vec3& v1, float k) {
  return Vec3(v1[0] + k, v1[1] + k, v1[2] + k);
}
inline Vec3 operator+(float k, const Vec3& v2) { return v2 + k; }

inline Vec3 operator-(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
}
inline Vec3 operator-(const Vec3& v1, float k) {
  return Vec3(v1[0] - k, v1[1] - k, v1[2] - k);
}
inline Vec3 operator-(float k, const Vec3& v2) {
  return Vec3(k - v2[0], k - v2[1], k - v2[2]);
}

inline Vec3 operator*(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2]);
}
inline Vec3 operator*(const Vec3& v1, float k) {
  return Vec3(v1[0] * k, v1[1] * k, v1[2] * k);
}
inline Vec3 operator*(float k, const Vec3& v2) { return v2 * k; }

inline Vec3 operator/(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1[0] / v2[0], v1[1] / v2[1], v1[2] / v2[2]);
}
inline Vec3 operator/(const Vec3& v1, float k) {
  return Vec3(v1[0] / k, v1[1] / k, v1[2] / k);
}
inline Vec3 operator/(float k, const Vec3& v2) {
  return Vec3(k / v2[0], k / v2[1], k / v2[2]);
}

inline float dot(const Vec3& v1, const Vec3& v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

inline Vec3 cross(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2],
              v1[0] * v2[1] - v1[1] * v2[0]);
}

inline float length(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline float length2(const Vec3& v) { return dot(v, v); }
inline Vec3 normalize(const Vec3& v) { return v / length(v); }

inline Vec3 worldToLocal(const Vec3& v, const Vec3& lx, const Vec3& ly,
                         const Vec3& lz) {
  return Vec3(dot(v, lx), dot(v, ly), dot(v, lz));
}

inline Vec3 localToWorld(const Vec3& v, const Vec3& lx, const Vec3& ly,
                         const Vec3& lz) {
  Vec3 ret;
  for (int i = 0; i < 3; ++i) {
    ret[i] = v[0] * lx[i] + v[1] * ly[i] + v[2] * lz[i];
  }
  return ret;
}

inline Vec3 sphericalToCartesian(float phi, float theta) {
  return Vec3(std::cos(phi) * std::sin(theta), std::cos(theta),
              std::sin(phi) * std::sin(theta));
}

struct Ray {
  Vec3 origin;
  Vec3 direction;
  static constexpr float tmin = 1e-5;
  mutable float tmax = std::numeric_limits<float>::max();

  Ray() {}
  Ray(const Vec3& origin, const Vec3& direction)
      : origin(origin), direction(direction) {}

  Vec3 operator()(float t) const { return origin + t * direction; }
};

struct SurfaceInfo {
  Vec3 position;
  Vec3 normal;
  Vec3 dpdu;
  Vec3 dpdv;
};

// forward declaration
class Primitive;

struct IntersectInfo {
  float t;
  SurfaceInfo surfaceInfo;
  const Primitive* hitPrimitive;
};

#endif