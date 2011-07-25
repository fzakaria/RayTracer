#include "mesh.hpp"
#include <iostream>
#include <limits.h>

#define BOUNDING_DEBUG 0

Mesh::Mesh(const std::vector<Point3D>& verts,
           const std::vector< std::vector<int> >& faces)
  : Primitive(), m_verts(verts),
    m_faces(faces)
{
	//Step 1: calculate the normal and d of each plane
	for (int i = 0 ; i < faces.size() ; ++i)
	{
		Face currentFace = m_faces[i];
		if (currentFace.size() < 3) { std::cerr << "ERROR: Face defined with less than 3 vertices" << std::endl; }
		Vector3D V0V1 = verts[currentFace[1]]-verts[currentFace[0]];
		Vector3D V0V2 = verts[currentFace[2]]-verts[currentFace[0]];
		
		Vector3D normal = V0V1.cross(V0V2);
		normal.normalize();
		
		m_normals.push_back(normal); //save the normal at same index as it's corresponding faec
		
		Vector3D V0 = verts[currentFace[0]] - Point3D(0,0,0);
		double d = (-1*V0).dot(normal);
		m_d.push_back(d);
	}
			
	//bounding box
	//step 1: find minX, minY, minZ and maxX, maxY, maxZ
	double minX = INT_MAX;
	double minY = INT_MAX;
	double minZ = INT_MAX;
	
	double maxX = INT_MIN;
	double maxY = INT_MIN;
	double maxZ = INT_MIN	;
	for (int i = 0 ; i < m_verts.size() ; ++i)
	{
		Point3D currVertex = m_verts[i];
		
		if (currVertex[0] < minX) { minX = currVertex[0]; }
		if (currVertex[0] > maxX) { maxX = currVertex[0]; }
		
		if (currVertex[1] < minY) { minY = currVertex[1]; }
		if (currVertex[1] > maxY) { maxY = currVertex[1]; }
		
		if (currVertex[2] < minZ) { minZ = currVertex[2]; }
		if (currVertex[2] > maxZ) { maxZ = currVertex[2]; }		
	}
	Point3D maxP(maxX, maxY, maxZ);
	Point3D minP(minX, minY, minZ);
	Point3D center((maxX+minX)/2, (maxY+minY)/2 , (maxZ+minZ)/2);
	Vector3D radius = maxP - minP;

	m_boundingSphere = new NonhierSphere(center, radius.length()/2 );
	
}

Mesh::~Mesh()
{
	delete m_boundingSphere;
}

bool Mesh::IntersectRay(Ray &r,Intersection &isect)
{
	Intersection copyIsect = isect;
	bool boundingHit = m_boundingSphere->IntersectRay(r, copyIsect);

	if (!boundingHit) { return false; } //bounding test failed
	
	if (BOUNDING_DEBUG) { isect = copyIsect; return true; } //if option set let us draw the bounding sphere

	bool hit = false;
	for (int i = 0 ; i < m_faces.size() ; ++i)
	{
		Face currentFace = m_faces[i];
		double ND = m_normals[i].dot(r.getDirection());
		if (ND == 0 ) { continue; } //reject intersection, it's parallel
		
		double t = ((m_verts[currentFace[0]] - r.getOrigin()).dot(m_normals[i]))/ND;
		
		if (t <= 0) { continue; } //reject intersection, intersection behind origin
		
		if (t > isect.t) { continue; } //reject as closer intersection found already
		
		Point3D p = r.getOrigin() + (t*r.getDirection() );
		
		//now we want to iterate over the triangles of each face
		//std::cout << "face-t calculated: " << t << std::endl;
	
		Point3D vertex0 = m_verts[currentFace[0]];
		for (int j = 2; j < currentFace.size(); ++j)
		{
			Point3D vertex1 = m_verts[currentFace[j-1]];
			Point3D vertex2 = m_verts[currentFace[j]];
			Triangle newTriangle(vertex0, vertex1, vertex2);
			
			bool trihit = intersectTriangle(newTriangle, r);
			if (trihit)
			{	
				//std::cout << "Triangle hit" << std::endl;
				isect.normal = m_normals[i];
				isect.intersect = p;
				isect.t = t;
				hit = true;
				break;
			}//if
		}//for
		
		
	
	}
	return hit;
}

bool Mesh::intersectTriangle(Triangle & triangle, Ray & r)
{
	//compute gamma
	double a = triangle.v0[0] - triangle.v1[0];
	double b = triangle.v0[1] - triangle.v1[1];
	double c = triangle.v0[2] - triangle.v1[2];
	double d = triangle.v0[0] - triangle.v2[0];
	double e = triangle.v0[1] - triangle.v2[1];
	double f = triangle.v0[2] - triangle.v2[2];
	double g = r.getDirection()[0];
	double h = r.getDirection()[1];
	double i = r.getDirection()[2];
	double j = triangle.v0[0] - r.getOrigin()[0];
	double k = triangle.v0[1] - r.getOrigin()[1];
	double l = triangle.v0[2] - r.getOrigin()[2];
	
	double M  = a*( (e*i) - (h*f) ) + (b*( (g*f) - (d*i) )) + (c*( (d*h) - (e*g)));
	
	
	double t = -1*((f*((a*k) - (j*b) )) + (e*((j*c)- (a*l))) + (d*((b*l) - (k*c))));
	t = t / M;
	
	//std::cout << "\ttriangle-t calculated: " << t << std::endl;
			
	//bounding box
	//step 1: find minX, minY, minZ and maxX, maxY, maxZ
	double gamma = i*( (a*k) - (j*b) ) + ( h* ((j*c) - (a*l)) ) + (g*((b*l) - (k*c)));
	gamma = gamma / M;
	
	if (gamma < 0 || gamma > 1) { return false; }
	
	double beta = j*((e*i) - (h*f)) + (k*((g*f) - (d*i))) + (l*((d*h)-(e*g)));
	beta = beta / M;
	if (beta < 0 || beta > (1-gamma) ) { return false; }

	return true;
}

std::ostream& operator<<(std::ostream& out, const Mesh& mesh)
{
  std::cerr << "mesh({";
  for (std::vector<Point3D>::const_iterator I = mesh.m_verts.begin(); I != mesh.m_verts.end(); ++I) {
    if (I != mesh.m_verts.begin()) std::cerr << ",\n      ";
    std::cerr << *I;
  }
  std::cerr << "},\n\n     {";
  
  for (std::vector<Mesh::Face>::const_iterator I = mesh.m_faces.begin(); I != mesh.m_faces.end(); ++I) {
    if (I != mesh.m_faces.begin()) std::cerr << ",\n      ";
    std::cerr << "[";
    for (Mesh::Face::const_iterator J = I->begin(); J != I->end(); ++J) {
      if (J != I->begin()) std::cerr << ", ";
      std::cerr << *J;
    }
    std::cerr << "]";
  }
  std::cerr << "});" << std::endl;
  return out;
}
