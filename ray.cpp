#include "ray.hpp"
#include "scene.hpp"

//IOR starts at 1.0 for air

Ray::Ray() : m_Origin(Point3D(0,0,0)), m_Direction(Vector3D(0,0,0)), m_bounce(0)
{
}

Ray::Ray(const Point3D& a_Origin,const Vector3D& a_Dir ) : m_Origin(a_Origin), m_Direction(a_Dir), m_bounce(0)
{
}
