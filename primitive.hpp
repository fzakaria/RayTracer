#ifndef CS488_PRIMITIVE_HPP
#define CS488_PRIMITIVE_HPP
#include "ray.hpp"
#include "algebra.hpp"
class NonhierSphere;
class NonhierBox;
class Image;

class Primitive {
public:
  Primitive();
  virtual ~Primitive();
  virtual bool IntersectRay(Ray &r,Intersection &isect);
  virtual void apply_texture(Image * texture) { m_texture = texture; }
  virtual volatile Image* getTexture() { return m_texture; }
  virtual Colour getTextureColour(const Intersection & isect);
  virtual void apply_bump(Image * bump) { m_bump = bump; }
  
  static Colour getImageColour(volatile Image * image, const Intersection & isect);
  
protected:
  volatile Image * m_texture;
  volatile Image * m_bump;
};


class Sphere : public Primitive {
public:
  Sphere();
  virtual ~Sphere();
  virtual bool IntersectRay(Ray &r,Intersection &isect);
  virtual void apply_texture(Image * texture);
  virtual void apply_bump(Image * bump);
private:
  NonhierSphere * m_sphere;
};

class Cube : public Primitive {
public:
  Cube();
  virtual ~Cube();
  virtual bool IntersectRay(Ray &r,Intersection &isect);
  virtual void apply_texture(Image * texture);
  virtual void apply_bump(Image * bump);
private:
  NonhierBox * m_box;
};

class Cone: public Primitive {
public:
 Cone(double radius, double height);
 virtual ~Cone();
 virtual bool IntersectRay(Ray &r,Intersection &isect);
private:
 double m_radius;
 double m_height;
};


class Plane: public Primitive {
public:
 Plane(double width, double height);
 virtual ~Plane();
 virtual bool IntersectRay(Ray &r,Intersection &isect);
private:
 double m_width;
 double m_height;
};


class Cylinder: public Primitive {
public:
 Cylinder(double radius, double height);
 virtual ~Cylinder();
 virtual bool IntersectRay(Ray &r,Intersection &isect);
private:
 double m_radius;
 double m_height;
};

class Torus: public Primitive {
public:
 Torus(double innerRadius, double outerRadius);
 virtual ~Torus();
 virtual bool IntersectRay(Ray &r,Intersection &isect);
private:
 double m_innerRadius;
 double m_outerRadius;
};

class NonhierSphere : public Primitive {
public:
  NonhierSphere(const Point3D& pos, double radius)
    : m_pos(pos), m_radius(radius)
  {
  }
  virtual ~NonhierSphere();
  virtual bool IntersectRay(Ray &r,Intersection &isect);
  
private:
  Point3D m_pos;
  double m_radius;
};

//NOTE: You want to do the inverse transfrmation of the Transformation to put the ray into object space
class NonhierBox : public Primitive {
public:
  NonhierBox(const Point3D& pos, double size);
  
  virtual ~NonhierBox();
  virtual bool IntersectRay(Ray &r,Intersection &isect);

private:
  Point3D m_pos;
  double m_size;
};

#endif
