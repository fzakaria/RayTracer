#include "algebra2.hpp"
#include <math.h>

// Return a matrix to represent a counterclockwise rotation of "angle"
// degrees around the axis "axis", where "axis" is one of the
// characters 'x', 'y', or 'z'.
#define PI 3.14159265

Matrix4x4 rotation(double angle, char axis)
{
  Matrix4x4 r;
  switch (axis)
  {
  case 'X':
  case 'x':
  	r[1][1] = cos(angle*PI/180);
  	r[1][2] = -1*sin(angle*PI/180);
  	r[2][1] = sin(angle*PI/180);
  	r[2][2] = cos(angle*PI/180);
  	break;
  case 'Y':
  case 'y':
    	r[0][0] = cos(angle*PI/180);
  	r[0][2] = sin(angle*PI/180);
  	r[2][0] = -1*sin(angle*PI/180);
  	r[2][2] = cos(angle*PI/180);
  	break;
  case 'Z':
  case 'z':
  	r[0][0] = cos(angle*PI/180);
  	r[0][1] = -1*sin(angle*PI/180);
  	r[1][0] = sin(angle*PI/180);
  	r[1][1] = cos(angle*PI/180);
  	break;
  default:
  	break;
  }
  return r;
}

// Return a matrix to represent a displacement of the given vector.
Matrix4x4 translation(const Vector3D& displacement)
{
  Matrix4x4 t;
  t[0][3] = displacement[0];
  t[1][3] = displacement[1];
  t[2][3] = displacement[2];
  return t;
}

// Return a matrix to represent a nonuniform scale with the given factors.
Matrix4x4 scaling(const Vector3D& scale)
{
  Matrix4x4 s;
  s[0][0] = scale[0];
  s[1][1] = scale[1];
  s[2][2] = scale[2];
  return s;
}
