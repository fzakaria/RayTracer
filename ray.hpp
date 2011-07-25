#ifndef CS488_RAY_HPP
#define CS488_RAY_HPP

#include "algebra.hpp"

class GeometryNode;
class Material;

class Ray
{
public:
	Ray();
	Ray(const Point3D& a_Origin, const Vector3D& a_Dir );
	
	Vector3D getDirection() const { return m_Direction; }
	
	Point3D  getOrigin() const { return m_Origin; }
	
	int getBounce() const { return m_bounce; }
	
	void setBounce(int value) { m_bounce = value; }
	
private:
	Point3D m_Origin;
	Vector3D m_Direction;
	int m_bounce;
};

class Intersection {
public:
	Intersection()
	{
		t = 100000;
		closestNode = NULL;
	}
      double t;
      Vector3D normal;
      Point3D intersect;
      GeometryNode * closestNode;
      Material * material;
      Point2D uv;

};

#endif
