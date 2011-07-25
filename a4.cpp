#include "a4.hpp"
#include "image.hpp"
#include "ray.hpp"
#include <algorithm>
#include <math.h>
#include <time.h>
#include <pthread.h>

#define PI 3.14159265
#define EPSILON 0.001

#define MAX_REFLECTION_COUNT 5

int Renderer::NUM_THREADS = 16;

bool Renderer::CELLSHADING_ON = 0;

int Renderer::AA_VALUE = 1;

bool Renderer::DOF_ON = false;
int Renderer::dispersion = 2;
int Renderer::DoFComplexity = 1;
int Renderer::DOF = 1;

void Renderer::set_dof(int z_dist, int dispersion_rad, int num_rays)
{
	DOF_ON = true; //turn on DOF
	DOF = z_dist;
	dispersion = dispersion_rad >=0 ? dispersion_rad : 0;
	DoFComplexity = num_rays >=1 ? num_rays : 1;	
}


void a4_render(// What to render
               SceneNode* root,
               // Where to output the image
               const std::string& filename,
               // Image size
               int width, int height,
               // Viewing parameters
               const Point3D& eye, const Vector3D& view,
               const Vector3D& up, double fov,
               // Lighting parameters
               const Colour& ambient,
               const std::list<Light*>& lights
               )
{
  Renderer * renderer = new Renderer(root, filename, width, height, eye, view, up, fov, ambient, lights);

  std::cerr << "a4_render(" << root << ",\n     "
            << filename << ", " << width << ", " << height << ",\n     "
            << eye << ", " << view << ", " << up << ", " << fov << ",\n     "
            << ambient << ",\n     {";

  for (std::list<Light*>::const_iterator I = lights.begin(); I != lights.end(); ++I) {
    if (I != lights.begin()) std::cerr << ", ";
  }
  std::cerr << "});" << std::endl;
  
  std::cerr << " Option Settings: \n";
  std::cerr <<  "\t#ofThreads: "<<Renderer::NUM_THREADS<<"x\n";
  std::cerr <<  "\tAA: "<<Renderer::AA_VALUE<<"x\n";
  if (Renderer::DOF_ON)
  {
  	std::cerr <<  "\tDOF is ON , z-plane @: " << Renderer::DOF << " CircleOfConfision: " << Renderer::dispersion << " #ofRays: " << Renderer::DoFComplexity <<"\n";
  }
  else { std::cerr <<  "\tDOF is OFF\n";  } 
  std::cerr << "\tCellShading: "<< Renderer::CELLSHADING_ON << std::endl;
  

  pthread_t thread_id[Renderer::NUM_THREADS];
  Renderer::ThreadParams params[Renderer::NUM_THREADS];
  //START THREADS
  for (int i = 0 ; i < Renderer::NUM_THREADS ; ++i)
  {
  	params[i] = Renderer::ThreadParams(i, renderer);
  	pthread_create(&thread_id[i],NULL, &Renderer::Begin, (void *)(params+i));
  }
  //WAIT FOR THREADS TO FINISH
  for (int i = 0 ; i < Renderer::NUM_THREADS ; ++i)
  {
  	  pthread_join(thread_id[i], NULL);
  }
  //cleanup
  delete renderer;
}

               
Renderer::Renderer(// What to render
               SceneNode* root,
               // Where to output the image
               const std::string& filename,
               // Image size
               int width, int height,
               // Viewing parameters
               const Point3D& eye, const Vector3D& view,
               const Vector3D& up, double fov,
               // Lighting parameters
               const Colour& ambient,
               const std::list<Light*>& lights
               ) : m_root(root), m_filename(filename), m_width(width), m_height(height), m_eye(eye),
               m_view(view), m_up(up), m_fov(fov), m_ambient(ambient), m_lights(lights)
{ 
  srand ( time(NULL) );
  //precomupted thread constant values
   // For now, just make a sample image.
  m_img = new Image(m_width, m_height, 3);
  
  //we do this go around volatile  
  Image * grayImage = new Image();
  if (!grayImage->loadPng("./Textures/grayscale.png")) { std::cerr << "Failed reading grayscale texture!" << std::endl;}
  m_grayScale = grayImage;
  
  //calculate focal length 
  m_focalLength = m_height / (2 * tan(PI * m_fov/(4*180)));
  //calculate orthonormal basis for camera
  m_w = -1*m_view; m_w.normalize();
  m_u = (m_up.cross(m_w)); m_u.normalize();
  m_v = (-1*m_w).cross(m_u); //Note that the y-axis is calculated using the reversed z-axis, the image turned out to be   upside down without this adjustment. I love trial-and-error. :)
  
  
}

Renderer::~Renderer()
{
	//write to image for deleting it. Doing it here insures that all threads should be finished
	std::cout << "Saving Image..." << std::endl;
	m_img->savePng(m_filename);
	delete m_img;
	delete m_grayScale;
}

void * Renderer::Begin(void * param)
{
	ThreadParams * threadParams = (ThreadParams *) param;
	threadParams->renderer->render(threadParams->threadNum);
	return NULL;
}


void Renderer::render(int threadNum)
{
 
 if (threadNum >= Renderer::NUM_THREADS ) { std::cerr << "ThreadNum is Incorrect!" << std::endl; return; } //sanity check
 
 //for Multithreading, we simply stripe the image
 int y_min = (threadNum) * (m_height/NUM_THREADS);
 int y_max = (threadNum+1) * (m_height/NUM_THREADS);
 
 for (int y =y_min; y < y_max; y++) 
 {
    for (int x = 0; x < m_width; x++) 
    {
    //update backgroundColour
    Colour backgroundColour = Colour(0,0,(double)y/ m_height);
    Colour finalColour(0,0,0);
    
    for (int q = 0 ; q < AA_VALUE; ++q)
    {
    	for (int p = 0 ; p < AA_VALUE; ++p)
    	{
    
   		double finalXPixel = x;
    		double finalYPixel = y;
   		if (AA_VALUE > 1)
    		{
    			//this makes it stratified sampling!
   			double randX = (double)rand()/RAND_MAX;
    			double randY = (double)rand()/RAND_MAX;
    			finalXPixel = x + ((p+randX)/AA_VALUE);
    			finalYPixel = y + ((q+randY)/AA_VALUE);
   		 }
    
    		//calculate u and v (pixel to screen coords)
    		double viewX = (-1*m_width) + ((2*m_width)*(finalXPixel+0.5)/m_width);
    		double viewY = (-1*m_height) + ((2*m_height)*(finalYPixel+0.5)/m_height);
    
   		Vector3D direction = (-1*m_focalLength*m_w) + (viewX*m_u) +(viewY*m_v);
    		direction.normalize();
    		//setup base current Colour
    		Colour currentPixelColour(0,0,0);
    		Ray ray(m_eye, direction); //the ray we're going to fire
    
    		if (DOF_ON)
    		{ 
    			//lets calculate point of intersection for each ray
    		 	double timeIntersectDoF =  ( (Point3D(0,0,DOF) - ray.getOrigin()).dot(Vector3D(0,0,1)) ) / (ray.getDirection().dot(Vector3D(0,0,1)));
   		 	Point3D DoFIntersectPoint = ray.getOrigin() + (timeIntersectDoF*ray.getDirection());

    			//lets perform DoF   
    			for (int i = 0 ; i < DoFComplexity ; ++i)
    			{
    				 //lets calculate a jittered eye spot now
     				Vector3D offsetDoF;
    				offsetDoF[0] = ((double) dispersion / RAND_MAX)* (1.0f * rand());
     				offsetDoF[1] = ( (double) dispersion / RAND_MAX)* (1.0f * rand());
     				offsetDoF[2] = 0;

     				Point3D newRayOrigin = ray.getOrigin() + offsetDoF;
    
     				Ray rayDoF = Ray(newRayOrigin, DoFIntersectPoint - newRayOrigin);
    
     				Colour lightCalcColour = RayTrace(rayDoF, backgroundColour);
     				currentPixelColour = currentPixelColour + lightCalcColour;
    			}//for - DOF
    			currentPixelColour = (double)(1.0/(double)DoFComplexity) * currentPixelColour;
    			}
    		else
    		{
    			//normal rayTr40ace
    			currentPixelColour = RayTrace(ray, backgroundColour);
    		}
    
    		finalColour = finalColour + currentPixelColour;
    	}//for - p
    }//for 
    
    finalColour = finalColour*(1.0/(AA_VALUE*AA_VALUE));

   //no need to lock since x, y should never be be shared
    (*m_img)(x, y, 0) = finalColour.R();
    (*m_img)(x, y, 1) = finalColour.G();
    (*m_img)(x, y, 2) = finalColour.B();

    }
 }
 
 pthread_exit((void*) 0);

}

void Renderer::LightCalculation(const std::list<Light*>& lights, const Intersection & intersection, const Ray & ray, Colour& pixelColour)
{
    Vector3D v = -1 * ray.getDirection(); v.normalize(); //vector back towards the eye   
    PhongMaterial * material = (PhongMaterial *)intersection.material; //grab material
    if (material == NULL ) { std::cerr << "Material is NULL" << std::endl; return; } //sanity check
    
    Colour initialDiffuseColour = intersection.closestNode->getDiffuseColour(intersection); //will return either KD or texture value
    
    pixelColour = m_ambient*initialDiffuseColour; //initialize colour to be ambient value
    //lets iterate through all40 the lights and compute the light
    for (std::list<Light*>::const_iterator I = lights.begin(); I != lights.end(); ++I)
    	{
    	 Light* light = *I; //set current light
    	 
    	 Vector3D l = light->position - intersection.intersect; l.normalize();
    	 Vector3D h = v+l; h.normalize();
    	 Vector3D n = intersection.normal; //normal should be normalized already
    	 
    	 //Shadow-check
    	 Intersection shadowRecord;
    	 Point3D offsetIntersection = intersection.intersect + (EPSILON*l);
    	 Ray shadowRay(offsetIntersection, l); //the ray we're going to fire
    	 
    	 m_root->IntersectRay(shadowRay, shadowRecord);
    	 if (shadowRecord.t > EPSILON && shadowRecord.t < 100000) //we've hit something so we are shadow
    	 {
    	 	PhongMaterial * shadowMat = (PhongMaterial *)shadowRecord.material; //grab material
    	 	
    	 	if (!shadowMat->isDielectric())
    		{
  	 		continue;
  	 	}
    	 }
    	 //lets calculate attenuation
    	 double r = (light->position - intersection.intersect).length();
    	 double attenFactor = 1/(light->falloff[0] + (light->falloff[1]*r) + (light->falloff[2]*r*r));
    	 
    	 Colour lightColour = attenFactor*light->colour;
    	 
    	 double phongExponent = pow(std::max(0.0, n.dot(h)), material->getShininess() );
  	 
  	 if (CELLSHADING_ON)
  	 {
  	 	double v = std::max(0.0,n.dot(l));
 
    	 	Colour diffuseLightColour(0,0,0);
    	 
    		Intersection celIntersect;
    	 	celIntersect.uv = Point2D(0, v);
    	 	diffuseLightColour = Primitive::getImageColour(m_grayScale, celIntersect);

  	 	pixelColour = pixelColour + (initialDiffuseColour*diffuseLightColour);
  	 } 
  	 else
  	 {
  	 	pixelColour = pixelColour + (initialDiffuseColour*(lightColour)*std::max(0.0,n.dot(l)));
  	 }
  	 pixelColour = pixelColour + (material->getKS()*(lightColour)*phongExponent);		
 	
 	}
 	
}

bool refract(Vector3D direction , Vector3D normal, double nt , double nold, Vector3D & refractDirection)
{	
	double n = nold/nt;
	if (nold != 1.0)
	{
		if (nold != nt) { std::cout << "you have embedded objects, this is not handled properly yet!" << std::cout; }
		n = 1.0/nold;
		normal = -1*normal;
	}

	
	double cosI = -1*direction.dot(normal);
	double cosT2 = 1.0 - n*n*(1.0- cosI * cosI);
	if (cosT2 > 0.0)
	{
		refractDirection = (n*direction) + (n * cosI - sqrt(cosT2)) * normal;
		return true;
	}

	return false;
}


bool refract(Vector3D direction, Vector3D normal, double IOR, Vector3D & refract)
{
        float n=1.0f/IOR;

        if (direction.dot(normal) >= 0)
        {
        	n = IOR;
        	normal = -1*normal;
        }
    
        float cosTi= -1*direction.dot(normal);
        float sinTt2=(float)n*n*(1.0-cosTi*cosTi);
        
        if (sinTt2<=1.0)
        {
            refract=n*direction+(n*cosTi-sqrtf(1.0-sinTt2))*normal;
            refract.normalize();
	    return true;
        }
        return false; //we have Total Internal refraction
}//refract()

/*We need to pass initial background colour, if we'd like to set one up that is non static (dynamic based on x,y) for instance*/
Colour Renderer::RayTrace(Ray & ray, Colour & backgroundColour)
{
    Intersection its; //record for intersection
    Colour lightCalcColour = backgroundColour;
    m_root->IntersectRay(ray, its);
    if (its.closestNode != NULL)
    {
	 //lets calcualte the light
   	 LightCalculation(m_lights, its, ray, lightCalcColour);
   	 PhongMaterial * material = (PhongMaterial *)its.material; //grab material
   	 
   	 if (material->isDielectric() && (ray.getBounce() <  MAX_REFLECTION_COUNT))
   	 { 	
   	 	Vector3D refractDirection;
   	 	
	 	if (refract(ray.getDirection(), its.normal , material->getRefractiveIndex() , refractDirection))
	 	{
	 		Point3D refractStartPoint = its.intersect + (EPSILON*refractDirection);
	 		Ray refractRay(refractStartPoint, refractDirection);
	 		refractRay.setBounce(ray.getBounce() + 1);
	 		lightCalcColour =  0.1*lightCalcColour + RayTrace(refractRay, backgroundColour);
		}

   	 }

   	 Colour reflectiveColour = material->getKM();
   	 if (!(reflectiveColour == Colour(0,0,0)))
   	 { 
   	 	Vector3D reflectedDirection = ray.getDirection() - (2*(ray.getDirection().dot(its.normal))*its.normal);
   	 	Point3D reflectedStartPoint = its.intersect + (EPSILON*reflectedDirection);
   	 	Ray reflectedRay(reflectedStartPoint, reflectedDirection);
   		if (ray.getBounce() <  MAX_REFLECTION_COUNT)
   	 	{
   	 		reflectedRay.setBounce(ray.getBounce() + 1);
   	 		lightCalcColour = lightCalcColour + (reflectiveColour * RayTrace(reflectedRay, backgroundColour));
   	 	}
   	 }
    }
    
    return lightCalcColour;
}
