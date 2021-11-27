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

constexpr float RAY_EPS = 1e-5f;

inline float rad2deg(float rad) { return 180.0f * rad / PI; }
inline float deg2rad(float deg) { return deg / 180.0f * PI; }

template <typename T>
struct Vec2 {
  T v[2];

  Vec2() { v[0] = v[1] = 0; }
  Vec2(T x) { v[0] = v[1] = x; }
  Vec2(T x, T y) {
    v[0] = x;
    v[1] = y;
  }

  T operator[](int i) const { return v[i]; }
  T& operator[](int i) { return v[i]; }

  Vec2 operator-() const { return Vec2(-v[0], -v[1]); }

  Vec2& operator+=(const Vec2& v) {
    for (int i = 0; i < 2; ++i) {
      this->v[i] += v[i];
    }
    return *this;
  }
  Vec2& operator*=(const Vec2& v) {
    for (int i = 0; i < 2; ++i) {
      this->v[i] *= v[i];
    }
    return *this;
  }
  Vec2& operator/=(const Vec2& v) {
    for (int i = 0; i < 2; ++i) {
      this->v[i] /= v[i];
    }
    return *this;
  }
};

template <typename T>
inline Vec2<T> operator+(const Vec2<T>& v1, const Vec2<T>& v2) {
  return Vec2<T>(v1[0] + v2[0], v1[1] + v2[1]);
}
template <typename T>
inline Vec2<T> operator+(const Vec2<T>& v1, float k) {
  return Vec2<T>(v1[0] + k, v1[1] + k);
}
template <typename T>
inline Vec2<T> operator+(float k, const Vec2<T>& v2) {
  return v2 + k;
}

template <typename T>
inline Vec2<T> operator-(const Vec2<T>& v1, const Vec2<T>& v2) {
  return Vec2<T>(v1[0] - v2[0], v1[1] - v2[1]);
}
template <typename T>
inline Vec2<T> operator-(const Vec2<T>& v1, float k) {
  return Vec2<T>(v1[0] - k, v1[1] - k);
}
template <typename T>
inline Vec2<T> operator-(float k, const Vec2<T>& v2) {
  return Vec2<T>(k - v2[0], k - v2[1]);
}

template <typename T>
inline Vec2<T> operator*(const Vec2<T>& v1, const Vec2<T>& v2) {
  return Vec2<T>(v1[0] * v2[0], v1[1] * v2[1]);
}
template <typename T>
inline Vec2<T> operator*(const Vec2<T>& v1, float k) {
  return Vec2<T>(v1[0] * k, v1[1] * k);
}
template <typename T>
inline Vec2<T> operator*(float k, const Vec2<T>& v2) {
  return v2 * k;
}

template <typename T>
inline Vec2<T> operator/(const Vec2<T>& v1, const Vec2<T>& v2) {
  return Vec2<T>(v1[0] / v2[0], v1[1] / v2[1]);
}
template <typename T>
inline Vec2<T> operator/(const Vec2<T>& v1, float k) {
  return Vec2<T>(v1[0] / k, v1[1] / k);
}
template <typename T>
inline Vec2<T> operator/(float k, const Vec2<T>& v2) {
  return Vec2<T>(k / v2[0], k / v2[1]);
}

using Vec2f = Vec2<float>;

template <typename T>
struct Vec3 {
  T v[3];

  // implement Point
  static constexpr int dim = 3;

  Vec3() { v[0] = v[1] = v[2] = 0; }
  Vec3(T x) { v[0] = v[1] = v[2] = x; }
  Vec3(T x, T y, T z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }

  T operator[](int i) const { return v[i]; }
  T& operator[](int i) { return v[i]; }

  Vec3 operator-() const { return Vec3(-v[0], -v[1], -v[2]); }

  Vec3& operator+=(const Vec3& v) {
    for (int i = 0; i < 3; ++i) {
      this->v[i] += v[i];
    }
    return *this;
  }
  Vec3& operator*=(const Vec3& v) {
    for (int i = 0; i < 3; ++i) {
      this->v[i] *= v[i];
    }
    return *this;
  }
  Vec3& operator/=(const Vec3& v) {
    for (int i = 0; i < 3; ++i) {
      this->v[i] /= v[i];
    }
    return *this;
  }
};

template <typename T>
inline Vec3<T> operator+(const Vec3<T>& v1, const Vec3<T>& v2) {
  return Vec3<T>(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
}
template <typename T>
inline Vec3<T> operator+(const Vec3<T>& v1, float k) {
  return Vec3<T>(v1[0] + k, v1[1] + k, v1[2] + k);
}
template <typename T>
inline Vec3<T> operator+(float k, const Vec3<T>& v2) {
  return v2 + k;
}

template <typename T>
inline Vec3<T> operator-(const Vec3<T>& v1, const Vec3<T>& v2) {
  return Vec3<T>(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
}
template <typename T>
inline Vec3<T> operator-(const Vec3<T>& v1, float k) {
  return Vec3<T>(v1[0] - k, v1[1] - k, v1[2] - k);
}
template <typename T>
inline Vec3<T> operator-(float k, const Vec3<T>& v2) {
  return Vec3<T>(k - v2[0], k - v2[1], k - v2[2]);
}

template <typename T>
inline Vec3<T> operator*(const Vec3<T>& v1, const Vec3<T>& v2) {
  return Vec3<T>(v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2]);
}
template <typename T>
inline Vec3<T> operator*(const Vec3<T>& v1, float k) {
  return Vec3<T>(v1[0] * k, v1[1] * k, v1[2] * k);
}
template <typename T>
inline Vec3<T> operator*(float k, const Vec3<T>& v2) {
  return v2 * k;
}

template <typename T>
inline Vec3<T> operator/(const Vec3<T>& v1, const Vec3<T>& v2) {
  return Vec3<T>(v1[0] / v2[0], v1[1] / v2[1], v1[2] / v2[2]);
}
template <typename T>
inline Vec3<T> operator/(const Vec3<T>& v1, float k) {
  return Vec3<T>(v1[0] / k, v1[1] / k, v1[2] / k);
}
template <typename T>
inline Vec3<T> operator/(float k, const Vec3<T>& v2) {
  return Vec3<T>(k / v2[0], k / v2[1], k / v2[2]);
}

template <typename T>
inline T dot(const Vec3<T>& v1, const Vec3<T>& v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

template <typename T>
inline Vec3<T> cross(const Vec3<T>& v1, const Vec3<T>& v2) {
  return Vec3<T>(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2],
                 v1[0] * v2[1] - v1[1] * v2[0]);
}

using Vec3f = Vec3<float>;
using Vec3ui = Vec3<uint32_t>;

inline float length(const Vec3f& v) { return std::sqrt(dot(v, v)); }
inline float length2(const Vec3f& v) { return dot(v, v); }
inline Vec3f normalize(const Vec3f& v) { return v / length(v); }

inline void orthonormalBasis(const Vec3f& n, Vec3f& t, Vec3f& b) {
  if (std::abs(n[1]) < 0.9f) {
    t = normalize(cross(n, Vec3f(0, 1, 0)));
  } else {
    t = normalize(cross(n, Vec3f(0, 0, -1)));
  }
  b = normalize(cross(t, n));
}

// transform direction from world to local
inline Vec3f worldToLocal(const Vec3f& v, const Vec3f& lx, const Vec3f& ly,
                          const Vec3f& lz) {
  return Vec3f(dot(v, lx), dot(v, ly), dot(v, lz));
}

// transform direction from local to world
inline Vec3f localToWorld(const Vec3f& v, const Vec3f& lx, const Vec3f& ly,
                          const Vec3f& lz) {
  Vec3f ret;
  for (int i = 0; i < 3; ++i) {
    ret[i] = v[0] * lx[i] + v[1] * ly[i] + v[2] * lz[i];
  }
  return ret;
}

// compute cartesian coordinates from spherical coordinates
inline Vec3f sphericalToCartesian(float theta, float phi) {
  return Vec3f(std::cos(phi) * std::sin(theta), std::cos(theta),
               std::sin(phi) * std::sin(theta));
}

struct Ray {
  Vec3f origin;
  Vec3f direction;
  static constexpr float tmin = RAY_EPS;
  mutable float tmax = std::numeric_limits<float>::max();

  Ray() {}
  Ray(const Vec3f& origin, const Vec3f& direction)
      : origin(origin), direction(direction) {}

  Vec3f operator()(float t) const { return origin + t * direction; }
};

struct SurfaceInfo {
  Vec3f position;
  Vec3f geometricNormal;
  Vec3f shadingNormal;
  Vec3f dpdu;  // tangent vector
  Vec3f dpdv;  // bitangent vector
  Vec2f texcoords;
  Vec2f barycentric;
};

// forward declaration
class Primitive;

struct IntersectInfo {
  float t;  // distance to the hit point
  SurfaceInfo surfaceInfo;
  const Primitive* hitPrimitive;
};

#endif