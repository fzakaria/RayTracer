#ifndef CS488_MATERIAL_HPP
#define CS488_MATERIAL_HPP

#include "algebra.hpp"

class Material {
public:
  virtual ~Material();
  virtual void apply_gl() const = 0;

protected:
  Material()
  {
  }
};

class PhongMaterial : public Material {
public:
  PhongMaterial(const Colour& kd, const Colour& ks, double shininess);
  virtual ~PhongMaterial();

  virtual void apply_gl() const;
  
  virtual Colour getKD() { return m_kd; } 
  
  virtual Colour getKS() { return m_ks; }
  
  virtual Colour getKM() { return m_km; }
  
  virtual void setKM(Colour km) { m_km = km; }
  
  virtual double getShininess() { return m_shininess; }
  
  virtual bool isDielectric() { return (m_refractiveIndex != 1.00) ;}
  
  virtual void setRefractiveIndex(double index) { m_refractiveIndex = index; }
  
  virtual double getRefractiveIndex() { return m_refractiveIndex; }

private:
  Colour m_kd;
  Colour m_ks;
  Colour m_km;
  
  double m_refractiveIndex;

  double m_shininess;
};


#endif
