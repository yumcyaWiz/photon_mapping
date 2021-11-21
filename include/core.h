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
  Vec2(float x) { v[0] = v[1] = x; }
  Vec2(float x, float y) {
    v[0] = x;
    v[1] = y;
  }

  float operator[](int i) const { return v[i]; }
  float& operator[](int i) { return v[i]; }

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
using Vec2i = Vec2<int>;
using Vec2ui = Vec2<uint32_t>;

struct Vec3 {
  float v[3];

  // implement Point
  static constexpr int dim = 3;

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

inline void orthonormalBasis(const Vec3& n, Vec3& t, Vec3& b) {
  if (std::abs(n[1]) < 0.9f) {
    t = normalize(cross(n, Vec3(0, 1, 0)));
  } else {
    t = normalize(cross(n, Vec3(0, 0, -1)));
  }
  b = normalize(cross(t, n));
}

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

inline Vec3 sphericalToCartesian(float theta, float phi) {
  return Vec3(std::cos(phi) * std::sin(theta), std::cos(theta),
              std::sin(phi) * std::sin(theta));
}

struct Ray {
  Vec3 origin;
  Vec3 direction;
  static constexpr float tmin = RAY_EPS;
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
  Vec2f texcoords;
  Vec2f barycentric;
};

// forward declaration
class Material;

struct IntersectInfo {
  float t;
  uint32_t primID;
  SurfaceInfo surfaceInfo;
  const Material* material;
};

#endif