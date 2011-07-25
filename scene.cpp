#include "scene.hpp"
#include <iostream>
#include "algebra2.hpp"

SceneNode::SceneNode(const std::string& name)
  : m_name(name)
{
}

SceneNode::~SceneNode()
{
}

void SceneNode::rotate(char axis, double angle)
{
  //std::cerr << "Stub: Rotate " << m_name << " around " << axis << " by " << angle << std::endl;
  Matrix4x4 newRotation = rotation(angle, axis);
  //std::cerr << "Old rotation:\n " << m_rotation << " new rotation:\n " << newRotation * m_rotation << std::endl;
  
  set_transform(m_trans*newRotation);
}

void SceneNode::scale(const Vector3D& amount)
{
  //std::cerr << "Stub: Scale " << m_name << " by " << amount << std::endl;
  Matrix4x4 newScale = scaling(amount);
  //std::cerr << "Old scale:\n " << m_scaling << " newScale:\n " << newScale * m_scaling << std::endl;
  
  set_transform(m_trans*newScale);
}

void SceneNode::translate(const Vector3D& amount)
{
  //std::cerr << "Stub: Translate " << m_name << " by " << amount << std::endl;
  Matrix4x4 newTrans = translation(amount);
  //std::cerr << "Old translate:\n " << m_translation << " new translate:\n " << newTrans * m_translation << std::endl;
  
  set_transform(m_trans*newTrans);
}

bool SceneNode::IntersectRay(Ray &r,Intersection &isect)
{
  Point3D newP = m_invtrans * r.getOrigin();
  Vector3D newD = m_invtrans * r.getDirection();
  Ray newRay(newP, newD);
  bool hit = false;
  for (ChildList::const_iterator i = m_children.begin() ; i != m_children.end() ; ++i)
  {
  	hit = hit | (*i)->IntersectRay(newRay, isect);
  	
  }
  
   if (hit)
  {
  	Point3D newIntersectionP = m_trans*isect.intersect;
  	Vector3D newNormal = m_invtrans.transpose() * isect.normal;
  	isect.intersect = newIntersectionP;
  	isect.normal = newNormal; isect.normal.normalize();
  }
  
  return hit;
  
}

bool SceneNode::is_joint() const
{
  return false;
}

JointNode::JointNode(const std::string& name)
  : SceneNode(name)
{
}

JointNode::~JointNode()
{
}

bool JointNode::is_joint() const
{
  return true;
}

void JointNode::set_joint_x(double min, double init, double max)
{
  m_joint_x.min = min;
  m_joint_x.init = init;
  m_joint_x.max = max;
}

void JointNode::set_joint_y(double min, double init, double max)
{
  m_joint_y.min = min;
  m_joint_y.init = init;
  m_joint_y.max = max;
}

GeometryNode::GeometryNode(const std::string& name, Primitive* primitive)
  : SceneNode(name),
    m_primitive(primitive)
{
	m_material = NULL;
}

GeometryNode::~GeometryNode()
{
}

bool GeometryNode::IntersectRay(Ray &r,Intersection &isect)
{
  Point3D newP = m_invtrans * r.getOrigin();
  Vector3D newD = m_invtrans * r.getDirection();
  Ray newRay(newP, newD);
  
  
  if (m_primitive)
  {
  	bool hit = m_primitive->IntersectRay(newRay, isect);
  	if (hit)
  	{
  	  isect.closestNode = this;
  	  isect.material = m_material;
  	 
  	  
  	 Point3D newIntersectionP = m_trans*isect.intersect;
  	 Vector3D newNormal = m_invtrans.transpose() * isect.normal;

  	 isect.intersect = newIntersectionP;
 	 isect.normal = newNormal; isect.normal.normalize();
 	 
  	}
  	
  	return hit;
  }
  
  for (ChildList::const_iterator i = m_children.begin() ; i != m_children.end() ; ++i)
  {
  	std::cerr << "This Geometry node has a child! It probably shouldn't!" << std::endl;
  	(*i)->IntersectRay(newRay, isect);
  }
  
  std::cerr << "This geometry node has no primitive? " << std::endl;
  return false; //if here then it had no primitive
  
}

Colour GeometryNode::getDiffuseColour(const Intersection &isect)
{
	Colour returnColour(1,0,0);
	if (m_primitive->getTexture() != NULL)
	{
		returnColour = m_primitive->getTextureColour(isect);
	}
	else
	{
		PhongMaterial * material = (PhongMaterial*) m_material;
		if (material == NULL) { std::cerr << "Invalid cast to PhongMaterial!" << std::endl; return returnColour; }
		
		returnColour = material->getKD();
	}
	return returnColour; 
}
 
