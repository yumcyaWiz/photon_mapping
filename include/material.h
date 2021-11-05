#ifndef _MATERIAL_H
#define _MATERIAL_H
#include "core.h"

class Material {
  virtual Vec3 sampleBRDF(const Vec3& wo, Vec3& wi, float& pdf) const;
};

class Diffuse : public Material {
 private:
  Vec3 rho;

  Diffuse(const Vec3& rho) : rho(rho) {}
};

class Mirror : public Material {};

class Glass : public Material {};

#endif