#ifndef CS488_MESH_HPP
#define CS488_MESH_HPP

#include <vector>
#include <iosfwd>
#include "primitive.hpp"
#include "algebra.hpp"

// A polygonal mesh.
class Mesh : public Primitive {
public:
  struct Triangle
  {
  	Triangle(Point3D p1, Point3D p2, Point3D p3): v0(p1), v1(p2), v2(p3) {}
  	Point3D v0;
  	Point3D v1;
  	Point3D v2;
  };
  
  Mesh(const std::vector<Point3D>& verts,
       const std::vector< std::vector<int> >& faces);
       
  ~Mesh();

  virtual bool IntersectRay(Ray &r,Intersection &isect);

  typedef std::vector<int> Face;
  
  bool intersectTriangle(Triangle & triangle, Ray & ray);
  
private:
  std::vector<Point3D> m_verts;
  std::vector<Face> m_faces;
  
  std::vector<Vector3D> m_normals;
  std::vector<double> m_d;
  
  NonhierSphere * m_boundingSphere;

  friend std::ostream& operator<<(std::ostream& out, const Mesh& mesh);
};

#endif
