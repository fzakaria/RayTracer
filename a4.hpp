#ifndef CS488_A4_HPP
#define CS488_A4_HPP

#include <string>
#include "algebra.hpp"
#include "scene.hpp"
#include "light.hpp"
#include "ray.hpp"

class Image;

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
               );
                         
class Renderer
{
public:
//setup LUA bindings for setup variables
static int NUM_THREADS;
static bool CELLSHADING_ON;
static int AA_VALUE;
static bool DOF_ON;
static int dispersion;
static int DoFComplexity;
static int DOF;

static void set_dof(int z_dist, int dispersion_rad, int num_rays);


struct ThreadParams
{
	int threadNum;
	Renderer *renderer;
	ThreadParams(){}
	ThreadParams(int & num, Renderer * render): threadNum(num), renderer(render) {}
};

Renderer(// What to render
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
               );
virtual ~Renderer();

static void * Begin(void * renderer);

void render(int threadNum);

void LightCalculation(const std::list<Light*>& lights, const Intersection & intersection, const Ray & ray, Colour & pixelColour);

Colour RayTrace(Ray & ray, Colour & backgroundColour);


private:
SceneNode* m_root;
std::string m_filename;
double m_width, m_height;
Point3D m_eye;
Vector3D m_view, m_up;
double m_fov;
Colour m_ambient;
std::list<Light*> m_lights;

//information each ray will need
//no need to mutex lock as we don't change the values (for threads)
double m_focalLength;
//w, u , v are for calculating the eye matrix
Vector3D m_w;
Vector3D m_u;
Vector3D m_v;
volatile Image * m_img;
volatile Image * m_grayScale;
};

#endif
