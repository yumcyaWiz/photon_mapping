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

inline float rad2deg(float rad) { return 180.0f * rad / PI; }
inline float deg2rad(float deg) { return deg / 180.0f * PI; }

struct Vec3 {
  float x;
  float y;
  float z;

  Vec3() { x = y = z = 0; }
  Vec3(float v) { x = y = z = v; }
  Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
};

inline Vec3 operator+(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}
inline Vec3 operator+(const Vec3& v1, float k) {
  return Vec3(v1.x + k, v1.y + k, v1.z + k);
}
inline Vec3 operator+(float k, const Vec3& v2) { return v2 + k; }

inline Vec3 operator-(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}
inline Vec3 operator-(const Vec3& v1, float k) {
  return Vec3(v1.x - k, v1.y - k, v1.z - k);
}
inline Vec3 operator-(float k, const Vec3& v2) {
  return Vec3(k - v2.x, k - v2.y, k - v2.z);
}

inline Vec3 operator*(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}
inline Vec3 operator*(const Vec3& v1, float k) {
  return Vec3(v1.x * k, v1.y * k, v1.z * k);
}
inline Vec3 operator*(float k, const Vec3& v2) { return v2 * k; }

inline Vec3 operator/(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}
inline Vec3 operator/(const Vec3& v1, float k) {
  return Vec3(v1.x / k, v1.y / k, v1.z / k);
}
inline Vec3 operator/(float k, const Vec3& v2) {
  return Vec3(k / v2.x, k / v2.y, k / v2.z);
}

inline float dot(const Vec3& v1, const Vec3& v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline Vec3 cross(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
              v1.x * v2.y - v1.y * v2.x);
}

inline float length(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline float length2(const Vec3& v) { return dot(v, v); }
inline Vec3 normalize(const Vec3& v) { return v / length(v); }

struct Ray {
  Vec3 origin;
  Vec3 direction;
  static constexpr float tmin = 1e-3;
  mutable float tmax = std::numeric_limits<float>::max();

  Ray() {}
  Ray(const Vec3& origin, const Vec3& direction)
      : origin(origin), direction(direction) {}

  Vec3 operator()(float t) const { return origin + t * direction; }
};

// forward declaration
class Primitive;

struct IntersectInfo {
  float t;
  Vec3 hitPos;
  Vec3 hitNormal;
  const Primitive* hitPrimitive;
};

#endif