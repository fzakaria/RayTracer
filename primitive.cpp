#include "primitive.hpp"
#include "polyroots.hpp"
#include "image.hpp"
#include <pthread.h>

#define PI 3.14159265358979323846

Primitive::Primitive()
{
	m_texture = NULL;
	m_bump = NULL;
}

Primitive::~Primitive()
{
}

bool Primitive::IntersectRay(Ray &,Intersection &){ return false;}

//Realistic ray tracing By Peter Shirley, R. Keith Morley
Colour Primitive::getImageColour(volatile Image * image, const Intersection & isect)
{	
	float u = isect.uv[0] - int(isect.uv[0]);
	float v = isect.uv[1] - int(isect.uv[1]);
	
	u *= (image->width() -3 );
	v *= (image->height() - 3);
	
	int iu = (int) u;
	int iv = (int) v;
	
	float tu = u - iu;
	float tv = v - iv;
	
	double r = (*image)(iu, iv, 0);
	double g = (*image)(iu, iv, 1);
	double b = (*image)(iu, iv, 2);
	Colour iuiv(r, g, b);
	
	r = (*image)(iu+1, iv, 0);
	g = (*image)(iu+1, iv, 1);
	b = (*image)(iu+1, iv, 2);
	Colour iu1iv(r, g, b);
	
	r = (*image)(iu, iv+1, 0);
	g = (*image)(iu, iv+1, 1);
	b = (*image)(iu, iv+1, 2);
	Colour iuiv1(r, g, b);
	
	r = (*image)(iu+1, iv+1, 0);
	g = (*image)(iu+1, iv+1, 1);
	b = (*image)(iu+1, iv+1, 2);
	Colour iu1iv1(r, g, b);
	
	//bilinear filtering
	Colour returnColour = iuiv*(1-tu)*(1-tv) + iu1iv*tu*(1-tv) + iuiv1*(1-tu)*tv + iu1iv1*tu*tv;
	
	return returnColour;
}

Colour Primitive::getTextureColour(const Intersection &isect)
{
	//sanity check
	if (m_texture == NULL)
	{
		std::cerr << "Your texture is non-existant!" << std::endl;
		return Colour(0,0,0);
	}
	
	return getImageColour(m_texture, isect);
}


Sphere::Sphere()
{
	m_sphere = new NonhierSphere(Point3D(0,0,0), 1);
}

Sphere::~Sphere()
{
	delete m_sphere;
}

void Sphere::apply_texture(Image * texture)
{
	m_texture = texture;
	m_sphere->apply_texture(texture);
}

void Sphere::apply_bump(Image * bump)
{
	m_bump = bump;
	m_sphere->apply_bump(bump);
}


bool Sphere::IntersectRay(Ray & r ,Intersection & isect)
{
	return m_sphere->IntersectRay(r, isect);
}

Cube::Cube()
{
	m_box = new NonhierBox(Point3D(0,0,0), 1);
}

Cube::~Cube()
{
	delete m_box;
}

void Cube::apply_texture(Image * texture)
{
	m_box->apply_texture(texture);
}

void Cube::apply_bump(Image * texture)
{
	m_box->apply_bump(texture);
}

bool Cube::IntersectRay(Ray & r ,Intersection & isect)
{
	return m_box->IntersectRay(r, isect);
}

Cone::Cone(double radius, double height) : m_radius(radius), m_height(height)
{
}

Cone::~Cone()
{
}


bool Cone::IntersectRay(Ray & ray ,Intersection & isect)
{
	return false; //TO DO !
	bool hit = false;
	double px2 = ray.getOrigin()[0]*ray.getOrigin()[0];
	double py2 = ray.getOrigin()[1]*ray.getOrigin()[1];
	double pz2 = ray.getOrigin()[2]*ray.getOrigin()[2];
	double vx2 =  ray.getDirection()[0]*ray.getDirection()[0];
	double vy2 = ray.getDirection()[1]*ray.getDirection()[1];
	double vz2 = ray.getDirection()[2]*ray.getDirection()[2];
	double r2h2 = (m_radius*m_radius)/(m_height*m_height);
	
	double pxvx =  ray.getOrigin()[0]*ray.getDirection()[0];
	double pyvy =  ray.getOrigin()[1]*ray.getDirection()[1];
	double pzvz =  ray.getOrigin()[2]*ray.getDirection()[2];
	
	double A = vx2 + vy2 - vz2;
	
	double B = 2*(pxvx + pyvy - pzvz);
	
	double C = px2 + py2 - pz2;
	
	double disc = (B*B) - (4*A*C);
	if (disc < 0) { return false; } //no hit

	double roots[2];
	size_t numHit = quadraticRoots(A, B, C, roots);
	
		//try side face
	for (size_t i = 0 ; i < numHit; ++i)
	{
		if (roots[i] < isect.t && roots[i] > 0) //found something closer
		{
			Point3D p = ray.getOrigin() + (roots[i]*ray.getDirection());
			if (p[2] >= 0 && p[2] <= m_height)
			{
				isect.t = roots[i];
				hit = true;
				//calculate the normal
				//step 1: calculate point p of intersection
				isect.normal =  p - Point3D(0,0, p[2]);	isect.normal.normalize();		
				isect.intersect = p; //store the intersection
			}
		}
	}
	
	return hit;
	
	
}

Plane::Plane(double width, double height) : m_width(width) , m_height(height)
{
}

Plane::~Plane()
{
}

bool Plane::IntersectRay(Ray & ray ,Intersection & isect)
{
	bool hit = false;
	//top face first
	Vector3D normal(0,1,0);
	Point3D pointOnPlane(0,0, 0);
	double ND = normal.dot(ray.getDirection());
	if (ND != 0 )
	{
		double t = ((pointOnPlane - ray.getOrigin()).dot(normal))/ND;
		if (t < isect.t && t > 0)
		{
			Point3D p = ray.getOrigin() + (t*ray.getDirection());
			if (((p[0] <= m_width) && (p[2] <= m_height)) && ((p[0] >= 0) && (p[2] >= 0)))
			{
				isect.t = t;
				hit = true;
				//calculate the normal
				//step 1: calculate point p of intersection
				isect.normal =  normal;	
				isect.intersect = p; //store the intersection
				
				double u = ( p[0]) / (m_width);
			 	double v = ( p[2]) / (m_height);
			 	isect.uv = Point2D(u,v);
				
			}
		}
	}
	
	//bottom plane
	normal = Vector3D(0,-1,0);
	pointOnPlane = Point3D(0,0, 0);
	ND = normal.dot(ray.getDirection());
	if (ND != 0 )
	{
		double t = ((pointOnPlane - ray.getOrigin()).dot(normal))/ND;
		if (t < isect.t && t > 0)
		{
			Point3D p = ray.getOrigin() + (t*ray.getDirection());
			if (((p[0] <= m_width) && (p[2] <= m_height)) && ((p[0] >= 0) && (p[2] >= 0)))
			{
				isect.t = t;
				hit = true;
				//calculate the normal
				//step 1: calculate point p of intersection
				isect.normal =  normal;	
				isect.intersect = p; //store the intersection
				
				double u = ( p[0]) / (m_width);
			 	double v = ( p[2]) / (m_height);
			 	isect.uv = Point2D(u,v);
			
			}
		}
	}
	
	if (m_bump != NULL && hit ) //lets pertube the normal
	{
		Colour bumpColour = getImageColour(m_bump, isect);
		Vector3D bumpNormal(2.0 * (bumpColour.R() - 0.5), 2.0 * (bumpColour.G() - 0.5), 2.0 * (bumpColour.B() - 0.5));
		isect.normal = isect.normal + bumpNormal;			  
	}
	return hit;
}


Cylinder::Cylinder(double r, double h) : m_radius(r), m_height(h)
{
}

Cylinder::~Cylinder()
{
}

bool Cylinder::IntersectRay(Ray & ray ,Intersection & isect)
{
	bool hit = false;
	
	double px2 = ray.getOrigin()[0]*ray.getOrigin()[0];
	double pz2 = ray.getOrigin()[2]*ray.getOrigin()[2];
	double vx2 = ray.getDirection()[0]*ray.getDirection()[0];
	double vz2 = ray.getDirection()[2]*ray.getDirection()[2];
	double pxvx = ray.getOrigin()[0]*ray.getDirection()[0];
	double pzvz = ray.getOrigin()[2]*ray.getDirection()[2];
	
	double A = vx2 + vz2;
	
	double B = 2*(pxvx + pzvz);
	
	double C = px2+ pz2 - (m_radius*m_radius);
	
	double disc = (B*B) - (4*A*C);
	if (disc < 0) { return false; } //no hit

	double roots[2];
	size_t numHit = quadraticRoots(A, B, C, roots);
	
	//try side face
	for (size_t i = 0 ; i < numHit; ++i)
	{
		if (roots[i] < isect.t && roots[i] > 0) //found something closer
		{
			Point3D p = ray.getOrigin() + (roots[i]*ray.getDirection());
			if (p[1] >= 0 && p[1] <= m_height)
			{
				isect.t = roots[i];
				hit = true;
				//calculate the normal
				//step 1: calculate point p of intersection
				isect.normal =  p - Point3D(0,p[1], 0);	isect.normal.normalize();		
				isect.intersect = p; //store the intersection
				
		
				isect.uv = Point2D(isect.normal[0]/2 + 0.5, isect.intersect[1]/m_height);
				
			}
		}
	}
	
	//try top and bottom faces
	//top first
	Vector3D normal(0,1,0);
	Point3D pointOnPlane(0,m_height, 0);
	double ND = normal.dot(ray.getDirection());
	if (ND != 0 )
	{
		double t = ((pointOnPlane - ray.getOrigin()).dot(normal))/ND;
		if (t < isect.t && t > 0)
		{
			Point3D p = ray.getOrigin() + (t*ray.getDirection());
			if ((p[0]*p[0]) + (p[2]*p[2]) <= (m_radius*m_radius))
			{
				isect.t = t;
				hit = true;
				//calculate the normal
				//step 1: calculate point p of intersection
				isect.normal =  normal;	
				isect.intersect = p; //store the intersection
				
				double u = (p[0] + m_radius)/(2*m_radius);
				double v = (p[2] + m_radius)/(2*m_radius);
				
				isect.uv = Point2D(u,v);
				
			}
		}
	}
	
	//bottom plane
	normal = Vector3D(0,-1,0);
	pointOnPlane = Point3D(0,0, 0);
	ND = normal.dot(ray.getDirection());
	if (ND != 0 )
	{
		double t = ((pointOnPlane - ray.getOrigin()).dot(normal))/ND;
		if (t < isect.t && t > 0)
		{
			Point3D p = ray.getOrigin() + (t*ray.getDirection());
			if ((p[0]*p[0]) + (p[2]*p[2]) <= (m_radius*m_radius))
			{
				isect.t = t;
				hit = true;
				//calculate the normal
				//step 1: calculate point p of intersection
				isect.normal =  normal;	
				isect.intersect = p; //store the intersection
				
				double u = (p[0] + m_radius)/(2*m_radius);
				double v = (p[2] + m_radius)/(2*m_radius);
				
				isect.uv = Point2D(u,v);
			}
		}
	}
	
	if (m_bump != NULL && hit) //lets pertube the normal
	{
				  Colour bumpColour = getImageColour(m_bump, isect);
				  Vector3D bumpNormal(2.0 * (bumpColour.R() - 0.5), 2.0 * (bumpColour.G() - 0.5), 2.0 * (bumpColour.B() - 0.5));
				  isect.normal = isect.normal + bumpNormal;			  
	}

	return hit;

}

Torus::Torus(double inner, double outer) : m_innerRadius(inner), m_outerRadius(outer)
{
}

Torus::~Torus()
{
}

//www.emeyex.com/site/projects/raytorus.pdf
bool Torus::IntersectRay(Ray & ray ,Intersection & isect)
{
	bool hit = false;
	double outerSquared = m_outerRadius*m_outerRadius;
	double innerSquared = m_innerRadius*m_innerRadius;
	
	double alpha = ray.getDirection().dot(ray.getDirection() );
	Vector3D originAsVector = ray.getOrigin() - Point3D(0,0,0);
	double beta = 2* ( originAsVector.dot(ray.getDirection()) );
	double gamma = (originAsVector.dot(originAsVector)) - (innerSquared) - (outerSquared);
	
	double A = alpha*alpha;
	double B = 2*alpha*beta;
	double C = (beta*beta) + (2*alpha*gamma)+ (4*outerSquared*ray.getDirection()[2]*ray.getDirection()[2]);
	double D = (2*beta*gamma) + (8*outerSquared*ray.getOrigin()[2]*ray.getDirection()[2]);
	double E = (gamma*gamma) + (4*outerSquared*ray.getOrigin()[2]*ray.getOrigin()[2]) - (4*outerSquared*innerSquared);
	
	double roots[4];
	size_t numRoots = quarticRoots(B/A, C/A, D/A, E/A, roots);
	
	for (size_t i = 0 ; i < numRoots; ++i)
	{
		if (roots[i] < isect.t && roots[i] > 0) //found something closer
		{
			isect.t = roots[i];
			hit = true;
			Point3D p = ray.getOrigin() + (roots[i] * ray.getDirection() );
			isect.intersect = p; //store point of intersection
			
			Vector3D normal;
			double x2 = p[0] * p[0];
			double y2 = p[1] * p[1];
			double z2 = p[2] * p[2];
			
			normal[0] = 4*p[0]*(x2 + y2 + z2 - innerSquared - outerSquared);
			normal[1] = 4*p[1]*(x2 + y2 + z2 - innerSquared - outerSquared);
			normal[2] = (4*p[2]*(x2 + y2 + z2 - innerSquared - outerSquared)) + (8*outerSquared*p[2]);
			isect.normal = normal;
			
			// Determine its angle from the y-axis.
			double u = (1.0 - (atan2(p[1], p[0]) + PI) / (PI*2));

			double len = sqrt(x2 + y2);

			// Now rotate about the y-axis to get the point P into the x-z plane.
			double x = len - m_outerRadius;
			double v = (atan2(p[2], x) + PI) / (PI*2);
			
			isect.uv = Point2D(u, v);
			
			if (m_bump != NULL ) //lets pertube the normal
			{
				  Colour bumpColour = getImageColour(m_bump, isect);
				  Vector3D bumpNormal(2.0 * (bumpColour.R() - 0.5), 2.0 * (bumpColour.G() - 0.5), 2.0 * (bumpColour.B() - 0.5));
				  isect.normal = isect.normal + bumpNormal;			  
			}

		}
	}
	
	return hit;
}



NonhierSphere::~NonhierSphere()
{
}

bool NonhierSphere::IntersectRay(Ray &r,Intersection &isect)
{
	bool hit = false;
	double A = r.getDirection().dot(r.getDirection() );
	Vector3D eMinusC = r.getOrigin() - m_pos;
	double B = (2*r.getDirection()).dot(eMinusC);
	double C = eMinusC.dot(eMinusC) - (m_radius*m_radius);
	
	double disc = (B*B) - (4*A*C);
	if (disc < 0) { return false; } //no hit
	
	double roots[2];
	size_t numHit = quadraticRoots(A, B, C, roots);
	for (size_t i = 0 ; i < numHit; ++i)
	{
		if (roots[i] < isect.t && roots[i] > 0)//found something closer
		{
			isect.t = roots[i];
			hit = true;
			//calculate the normal
			//step 1: calculate point p of intersection
			Point3D p = r.getOrigin() + (roots[i]*r.getDirection());
			isect.normal =  (1/m_radius)*(p - m_pos);

			isect.intersect = p; //store the intersection
			
			
			//calculate uv coordinates
			float theta = acos(isect.normal[2]);
			float phi = atan2(isect.normal[1], isect.normal[0]);
			
			if (phi < 0.0f) phi += 2*PI;

			isect.uv = Point2D(phi/(2*PI), (PI- theta)/PI);
			
			if (m_bump != NULL ) //lets pertube the normal
			{
				  Colour bumpColour = getImageColour(m_bump, isect);
				  Vector3D bumpNormal(2.0 * (bumpColour.R() - 0.5), 2.0 * (bumpColour.G() - 0.5), 2.0 * (bumpColour.B() - 0.5));
				  isect.normal = isect.normal + bumpNormal;			  
			}
		}
	}
	
	
	return hit;
}

NonhierBox::NonhierBox(const Point3D& pos, double size)
    : m_pos(pos), m_size(size)
{
}

NonhierBox::~NonhierBox()
{
}


/*The intersection test below was from: 3D math primer for graphics and game development
 By Fletcher Dunn, Ian Parberry*/
//EXTRA NOTE: For some reason , I had to comment out the early quits. Too tired to really think whhat they were used for.
//also if we were inside the box I just returned false immediately
bool NonhierBox::IntersectRay(Ray &r,Intersection &isect)
{

	//checkfor point in box, trivial reject, and determine parameteric distance to each front face
	//Point3D min( m_pos[0] - (m_size/2), m_pos[1] - (m_size/2), m_pos[2] - (m_size/2));
	//Point3D max( m_pos[0] + (m_size/2), m_pos[1] + (m_size/2), m_pos[2] + (m_size/2));
	Point3D min(m_pos[0], m_pos[1], m_pos[2]);
	Point3D max( m_pos[0] + (m_size), m_pos[1] + (m_size), m_pos[2] + (m_size));
	bool inside = true;
	Point3D rayOrigin = r.getOrigin();
	Vector3D rayDirection = r.getDirection();
	Vector3D normal;
	float xt, xn;
	bool hit = false;
	
	if (rayOrigin[0] < min[0])
	{
		xt = min[0] - rayOrigin[0];
		if (xt > rayDirection[0])
		{
		 //std::cout << "early quit #1" << std::endl;
		 //return false;
		}
		xt /= rayDirection[0];
		inside = false;
		xn = -1.0f;
	}
	else if (rayOrigin[0] > max[0] ) 
	{
		xt = max[0] - rayOrigin[0];
		if (xt < rayDirection[0])
		{
		 //std::cout << "early quit #2" << std::endl;
		 //return false;
		}
		xt /= rayDirection[0];
		inside = false;
		xn = 1.0f;
	}
	else
	{
		xt = -1.0f;
	}
	
	float yt, yn;
	if (rayOrigin[1] < min[1])
	{
		yt = min[1] - rayOrigin[1];
		if (yt > rayDirection[1])
		{
		 //std::cout << "early quit #3" << std::endl;
		 //return false;
		}
		yt /= rayDirection[1];
		inside = false;
		yn = -1.0f;
	}
	else if (rayOrigin[1] > max[1] ) 
	{
		yt = max[1] - rayOrigin[1];		if (m_bump != NULL ) //lets pertube the normal
		{
			Colour bumpColour = getImageColour(m_bump, isect);
			Vector3D bumpNormal(2.0 * (bumpColour.R() - 0.5), 2.0 * (bumpColour.G() - 0.5), 2.0 * (bumpColour.B() - 0.5));
			isect.normal = isect.normal + bumpNormal;			  
		}
		if (yt < rayDirection[1])
		{
		 //std::cout << "early quit #4" << std::endl;
		 //return false;
		}
		yt /= rayDirection[1];
		inside = false;
		yn = 1.0f;
	}
	else
	{
		yt = -1.0f;
	}
	
	float zt, zn;
	if (rayOrigin[2] < min[2])
	{
		zt = min[2] - rayOrigin[2];
		if (zt > rayDirection[2])
		{
		 //std::cout << "early quit #5" << std::endl;
		 //return false;
		}
		zt /= rayDirection[2];
		inside = false;
		zn = -1.0f;
	}
	else if (rayOrigin[2] > max[2] ) 
	{
		zt = max[2] - rayOrigin[2];
		if (zt < rayDirection[2])
		{
		 //std::cout << "early quit #6" << std::endl;
		 //return false;
		}
		zt /= rayDirection[2];
		inside = false;
		zn = 1.0f;
	}
	else
	{
		zt = -1.0f;
	}
	
	if (inside)
	{
		//this could be wrong: check if stuff messes up later
		//std::cout << "Inside! " << std::endl;
		isect.normal = -1*r.getDirection();
		isect.normal.normalize();
		isect.intersect = r.getOrigin();
		isect.t = 0;
		return true;
	}
	
	//select frathest plane = this is the plane of intersection
	int which = 0;
	float t = xt;
	if (yt > t){
		which = 1; t = yt;
	}
	if (zt > t) {
		which = 2; t = zt;
	}
	double u, v;
	switch (which){
	
		case 0: //intersect yz plane
		{
			float y = rayOrigin[1] + rayDirection[1]*t;
			if (y < min[1] || y > max[1]) return false;
			float z = rayOrigin[2] + rayDirection[2]*t;
			if (z < min[2] || z > max[2]) return false;
			normal = Vector3D(xn, 0, 0);
			
			 u = ( z - min[2]) / (max[2] - min[2]);
			 v = ( y - min[1]) / (max[1] - min[1]);
		}break;
		case 1: //intersect xz plane
		{
			float x = rayOrigin[0] + rayDirection[0]*t;
			if (x < min[0] || x > max[0]) return false;
			float z = rayOrigin[2] + rayDirection[2]*t;
			if (z < min[2] || z > max[2]) return false;
			normal = Vector3D(0, yn, 0);
			
			 u = ( x - min[0]) / (max[0] - min[0]);
			 v = ( z - min[2]) / (max[2] - min[2]);
		}break;
		case 2: //intersect xy plane
		{
			float x = rayOrigin[0] + rayDirection[0]*t;
			if (x < min[0] || x > max[0]) return false;
			float y = rayOrigin[1] + rayDirection[1]*t;
			if (y < min[1] || y > max[1]) return false;
			normal = Vector3D(0, 0, zn);
			
			 u = ( x - min[0]) / (max[0] - min[0]);
			 v = ( y - min[1]) / (max[1] - min[1]);
		}break;
	}
	
	if ( t > 0 && t < isect.t)
	{
		isect.t = t;
		hit = true;
		isect.normal = normal;
		Point3D p = r.getOrigin() + (t*r.getDirection());
		isect.intersect = p;
		
		isect.uv = Point2D(u, v);
		
		if (m_bump != NULL ) //lets pertube the normal
		{
			Colour bumpColour = getImageColour(m_bump, isect);
			Vector3D bumpNormal(2.0 * (bumpColour.R() - 0.5), 2.0 * (bumpColour.G() - 0.5), 2.0 * (bumpColour.B() - 0.5));
			isect.normal = isect.normal + bumpNormal;			  
		}
	}

	return hit;		
}
