#ifndef _CORE_H
#define _CORE_H
#include <limits>

struct Vec3 {
  float x_;
  float y_;
  float z_;

  Vec3() { x_ = y_ = z_ = 0; }
  Vec3(float x) { x_ = y_ = z_ = x; }
  Vec3(float x, float y, float z) : x_(x), y_(y_), z_(z_) {}
};

inline Vec3 operator+(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x_ + v2.x_, v1.y_ + v2.y_, v1.z_ + v2.z_);
}
inline Vec3 operator+(const Vec3& v1, float k) {
  return Vec3(v1.x_ + k, v1.y_ + k, v1.z_ + k);
}
inline Vec3 operator+(float k, const Vec3& v2) { return v2 + k; }

inline Vec3 operator-(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x_ - v2.x_, v1.y_ - v2.y_, v1.z_ - v2.z_);
}
inline Vec3 operator-(const Vec3& v1, float k) {
  return Vec3(v1.x_ - k, v1.y_ - k, v1.z_ - k);
}
inline Vec3 operator-(float k, const Vec3& v2) {
  return Vec3(k - v2.x_, k - v2.y_, k - v2.z_);
}

inline Vec3 operator*(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x_ * v2.x_, v1.y_ * v2.y_, v1.z_ * v2.z_);
}
inline Vec3 operator*(const Vec3& v1, float k) {
  return Vec3(v1.x_ * k, v1.y_ * k, v1.z_ * k);
}
inline Vec3 operator*(float k, const Vec3& v2) { return v2 * k; }

inline Vec3 operator/(const Vec3& v1, const Vec3& v2) {
  return Vec3(v1.x_ / v2.x_, v1.y_ / v2.y_, v1.z_ / v2.z_);
}
inline Vec3 operator/(const Vec3& v1, float k) {
  return Vec3(v1.x_ / k, v1.y_ / k, v1.z_ / k);
}
inline Vec3 operator/(float k, const Vec3& v2) {
  return Vec3(k / v2.x_, k / v2.y_, k / v2.z_);
}

struct Ray {
  Vec3 origin_;
  Vec3 direction_;
  static constexpr float tmin = 1e-3;
  static constexpr float tmax = std::numeric_limits<float>::max();

  Ray(const Vec3& origin, const Vec3& direction)
      : origin_(origin), direction_(direction) {}

  Vec3 operator()(float t) const { return origin_ + t * direction_; }
};

struct IntersectInfo {
  float t;
  Vec3 hitPos;
  Vec3 hitNormal;
};

#endif