/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cylinder class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#ifndef H_CYLINDER
#define H_CYLINDER
#include <glm/glm.hpp>
#include "SceneObject.h"

/**
 * Defines a simple Cylinder located at 'center'
 * with the specified radius
 */
class Cylinder : public SceneObject
{

private:
    glm::vec3 centerBottom;
    float radius;
    float height;

public:
	Cylinder()
		: centerBottom(glm::vec3(0)), radius(1.0), height(1.0)  
	{
		color = glm::vec3(1);
	};

    Cylinder(glm::vec3 c, float r, float height, glm::vec3 col)
		: centerBottom(c), radius(r), height(height)
	{
		color = col;
	};

	float intersect(glm::vec3 posn, glm::vec3 dir);

	glm::vec3 normal(glm::vec3 p);

};

#endif //!H_CYLINDER