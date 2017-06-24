/*----------------------------------------------------------
* COSC363  Ray Tracer Assignment
*
*  The cylinder class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cylinder.h"
#include "Plane.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include "TextureBMP.h"
/**
* cylinder's intersection method.  The input is a ray (pos, dir). 
*/
float Cylinder::intersect(glm::vec3 posn, glm::vec3 dir)
{
	float xc = centerBottom.x;
	float zc = centerBottom.z;
	float x = posn.x;
	float z = posn.z;
	
	float a = dir.x * dir.x + dir.z * dir.z;
	float b = 2 * (dir.x * (x - xc) + dir.z * (z - zc));
    float c = ((x - xc) * (x - xc)) + ((z - zc) * (z - zc)) - (radius * radius);
    float delta = b * b - 4 * (a * c);
   
	if(fabs(delta) < 0.001) return -1.0; 
    if(delta < 0.0) return -1.0;

    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);
    float t3 = (height - posn.y + centerBottom.y) / (dir.y);
    float t4 = (centerBottom.y - posn.y) / (dir.y);

	
	bool valid1 = true;
	bool valid2 = true;
	bool valid3 = false;
	bool valid4 = false;
    
    //check height bounds
    if (posn.y + t1 * dir.y - centerBottom.y < 0|| posn.y + t1 * dir.y - centerBottom.y > height){
		valid1 = false;
	}
	
	if (posn.y + t2 * dir.y - centerBottom.y < 0 || posn.y + t2 * dir.y - centerBottom.y > height){
		valid2 = false;
	}

	
	if (posn.y + t1 * dir.y - centerBottom.y >= height && posn.y + t2 * dir.y - centerBottom.y <= height){
		valid3 = true;
		if(fabs(t3) < 0.001 ){
			valid3 = false;
		}
	}

	if (posn.y + t1 * dir.y - centerBottom.y <= 0 && posn.y + t2 * dir.y - centerBottom.y >= 0){
		valid4 = true;
		if(fabs(t4) < 0.001 ){
			valid4 = false;
		}
	}
	
	if (valid3 && !valid4) { return t3; }
    if (!valid3 && valid4) { return t4;}
    if (valid3 && valid4) {return std::min(t3, t4);}
    
    if(fabs(t1) < 0.001 )
    {
        valid1 = false;
    }
    if(fabs(t2) < 0.001 ) valid2 = false;
    
    if (valid1 && !valid2) { return t1; }
    if (!valid1 && valid2) { return t2;}
    if (valid1 && valid2) {return std::min(t1, t2);}
    
	return -1;
}

/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the sphere.
*/
glm::vec3 Cylinder::normal(glm::vec3 p)
{
	if (p.y == height + centerBottom.y){
		return glm::vec3(0, 1, 0);
	}
	if (p.y == centerBottom.y){
		return glm::vec3(0, -1, 0);
	}
    glm::vec3 n = glm::vec3((p.x - centerBottom.x) / radius, 0.0f, (p.z - centerBottom.z) / radius);

    return n;
}
