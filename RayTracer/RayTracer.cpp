/*========================================================================
* COSC 363  Computer Graphics (2017)
* Ray tracer 
* See Lab07.pdf for details.
*=========================================================================
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "SceneObject.h"
#include "Ray.h"
#include <GL/glut.h>
#include "math.h"
#include "Plane.h"
#include "TextureBMP.h"
#include "Cylinder.h"
using namespace std;

const float WIDTH = 20.0;  
const float HEIGHT = 20.0;
const float EDIST = 40.0;
const int NUMDIV = 600;
const int MAX_STEPS = 5;
const float XMIN = -WIDTH * 0.5;
const float XMAX =  WIDTH * 0.5;
const float YMIN = -HEIGHT * 0.5;
const float YMAX =  HEIGHT * 0.5;

vector<SceneObject*> sceneObjects;  //A global list containing pointers to objects in the scene

TextureBMP texture;
TextureBMP wallTexture;
TextureBMP happyTexture;

//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
	glm::vec3 backgroundCol(0);
	glm::vec3 light(50, 40, -3);
	glm::vec3 light2(-30, 40, -3);
	float ambientTerm = 0.2;

    ray.closestPt(sceneObjects);		//Compute the closest point of intersetion of objects with the ray

    if(ray.xindex == -1) return backgroundCol;      //If there is no intersection return background colour
    
    glm::vec3 col = sceneObjects[ray.xindex]->getColor(); //object's colour;
    
    glm::vec3 normalVector = sceneObjects[ray.xindex]->normal(ray.xpt);
    glm::vec3 lightVector = light -ray.xpt;
    glm::vec3 lightVector2 = light2 - ray.xpt;
    glm::vec3 lightVectorNormal = normalize(lightVector);
    glm::vec3 lightVectorNormal2 = normalize(lightVector2);
    glm::vec3 reflVector = glm::reflect(-lightVectorNormal, normalVector);
    glm::vec3 reflVector2 = glm::reflect(-lightVectorNormal2, normalVector);
    
    if (ray.xindex == 3){
		float a1 = -20;
		float a2 = 20;
		float b1 = -40;
		float b2 = -200;
		float texcoords = (ray.xpt.x - a1) / (a2 - a1);
		float texcoordt = (ray.xpt.z - b1) / (b2 - b1);
		col = texture.getColorAt(texcoords, texcoordt);
	}
	
	if (ray.xindex == 2){
		// float texcoords = 0.5 + atan2(normalVector.z, normalVector.x) / (2*M_PI);
		float texcoords = 0.5 - asin(normalVector.x) / (M_PI);
		float texcoordt = 0.5 - asin(normalVector.y) / (M_PI);
		col = happyTexture.getColorAt(texcoords, texcoordt);	
	}
	
	 if (ray.xindex == 5){
		if (sin(ray.xpt.x * 3) > 0 && cos(ray.xpt.y * 3) < 0){
			col = glm::vec3(1, 0.5, 0);
		} else {
			col = glm::vec3(sin(ray.xpt.x), sin(ray.xpt.y), sin(ray.xpt.z));
		}
		
	}

    float lDotn = glm::dot(normalVector, lightVectorNormal);
    float rDotv = glm::dot(reflVector, -ray.dir);
    float lightDist = glm::length(light - ray.xpt);
    float lDotn2 = glm::dot(normalVector, lightVectorNormal2);
    float rDotv2 = glm::dot(reflVector2, -ray.dir);
    float lightDist2 = glm::length(light2 - ray.xpt);
    
    glm::vec3 specCol(1, 1, 1);
    glm::vec3 colorSum(0);
    float specTerm;
    float specTerm2;
    if (rDotv >= 0){
		specTerm = pow(rDotv, 10);
	} else{
		specTerm = 0;
	}
	
	if (rDotv2 >= 0){
		specTerm2 = pow(rDotv2, 5);
	} else{
		specTerm2 = 0;
	}

    Ray shadow(ray.xpt, lightVectorNormal);
	shadow.closestPt(sceneObjects);
	Ray shadow2(ray.xpt, lightVectorNormal2);
	shadow2.closestPt(sceneObjects);
	
	
	bool firstShadow = lDotn <= 0 || (shadow.xindex > -1 && shadow.xdist < lightDist);
	bool secondShadow = lDotn2 <= 0 || (shadow2.xindex > -1 && shadow2.xdist < lightDist2);

	if ((firstShadow && !secondShadow)){
		if (shadow.xindex == 6){
			colorSum = col * 0.3f * sceneObjects[shadow.xindex]->getColor();
		} else if (shadow.xindex == 1){
			colorSum = col * 0.4f;
		}
		else {
			colorSum = col * ambientTerm;
		}
		colorSum += lDotn2 * col;
	}else if ((!firstShadow && secondShadow)){
		if (shadow2.xindex == 6){
			colorSum = col * 0.4f * sceneObjects[shadow2.xindex]->getColor();
		}  else if (shadow2.xindex == 1){
			colorSum = col * 0.4f;
		}else {
			colorSum = col * ambientTerm;
		}
		colorSum += lDotn * col;
	} else if (!firstShadow && !secondShadow){
		colorSum = col * ambientTerm;
		colorSum += (lDotn * col + specCol * specTerm) + (lDotn2 * col + specCol * specTerm2);
	} else if (firstShadow && secondShadow){
		if (shadow.xindex == 6 || shadow2.xindex == 6){
			if (shadow.xindex > -1){
				colorSum = col * 0.4f * sceneObjects[shadow.xindex]->getColor();
				if (shadow2.xindex > -1){
					colorSum *= sceneObjects[shadow2.xindex]->getColor();
				}
			} else {
				colorSum = col * 0.4f * sceneObjects[shadow2.xindex]->getColor();
			}
		if (shadow.xindex == 1 || shadow2.xindex == 1){
			colorSum *= col * 0.4f;
		}
		} else {
			colorSum = col * ambientTerm;
		}
	}
	
	
	if(ray.xindex == 0 && step < MAX_STEPS)
	{
		glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
		Ray reflectedRay(ray.xpt, reflectedDir);
		glm::vec3 reflectedCol = trace(reflectedRay, step+1); //Recursion!
		colorSum = colorSum + (0.8f*reflectedCol);  
	}
	
	if(ray.xindex == 1 && step < MAX_STEPS)
	{
		float eta = 0.98f;
		glm::vec3 refactedDir = glm::refract(ray.dir, normalVector, eta);
		Ray refractedRay(ray.xpt, refactedDir);
		refractedRay.closestPt(sceneObjects);
		if (refractedRay.xindex > 0){
			glm::vec3 mNormalVector = sceneObjects[refractedRay.xindex]->normal(refractedRay.xpt);
			glm::vec3 mRefractedDirection = glm::refract(refactedDir, -mNormalVector, 1.0f / eta);
			Ray mrefractedRay(refractedRay.xpt, mRefractedDirection);
			glm::vec3 refractedCol = trace(mrefractedRay, step+1);
			colorSum = colorSum + (refractedCol);
		}
		
	}
	
	if(ray.xindex == 6 && step < MAX_STEPS)
	{
		float eta = 1.0f;
		glm::vec3 refactedDir = glm::refract(ray.dir, normalVector, eta);
		Ray refractedRay(ray.xpt, refactedDir);
		refractedRay.closestPt(sceneObjects);
		if (refractedRay.xindex > 0){
			glm::vec3 mNormalVector = sceneObjects[refractedRay.xindex]->normal(refractedRay.xpt);
			glm::vec3 mRefractedDirection = glm::refract(refactedDir, -mNormalVector, 1.0f / eta);
			Ray mrefractedRay(refractedRay.xpt, mRefractedDirection);
			glm::vec3 refractedCol = trace(mrefractedRay, step+1);
			colorSum = colorSum + (refractedCol);
		}
		
	}

	
	return colorSum;
}

bool compareVec3(glm::vec3 col1, glm::vec3 col2){
	float x = abs(col1.x - col2.x);
	float y = abs(col1.y - col2.y);
	float z = abs(col1.z - col2.z);
	float threshold = 0.25f;
	return x < threshold && y < threshold && z < threshold;
}

glm::vec3 antiAliasing(float xp, float yp, float cellX, float cellY, int step){
	glm::vec3 cols[4];
	glm::vec3 eye(0., 0., 0.);  //The eye position (source of primary rays) is the origin
	
	glm::vec3 dir1(xp+0.25*cellX, yp+0.25*cellY, -EDIST);	//direction of the left upper ray
	glm::vec3 dir2(xp+0.25*cellX, yp+0.75*cellY, -EDIST);	//direction of the right upper ray
	glm::vec3 dir3(xp+0.75*cellX, yp+0.25*cellY, -EDIST);	//direction of the left lower ray
	glm::vec3 dir4(xp+0.75*cellX, yp+0.75*cellY, -EDIST);	//direction of the right lower ray

	Ray ray1 = Ray(eye, dir1);		//Create a ray originating from the camera in the direction 'dir1'
	ray1.normalize();				//Normalize the direction of the ray to a unit vector
	glm::vec3 col = trace (ray1, 1); //Trace the primary ray and get the colour value
	cols[0] = col;

	Ray ray2 = Ray(eye, dir2);		
	ray2.normalize();				
	glm::vec3 col2 = trace (ray2, 1); 
	cols[1] = col2;

	Ray ray3 = Ray(eye, dir3);		
	ray3.normalize();		
	glm::vec3 col3 = trace (ray3, 1); 
	cols[2] = col3;
	
	Ray ray4 = Ray(eye, dir4);
	ray4.normalize();		
	glm::vec3 col4 = trace (ray4, 1); 	
	cols[3] = col4;
	glm::vec3 colTot = (col + col2 + col3 + col4) / 4.0f;
	step += 1;
	if (step < 10){
		if (!(compareVec3(col, colTot))){
			col = antiAliasing(xp, yp, cellX / 2.0f, cellY / 2.0f, step);
		}
		if (!(compareVec3(col2, colTot))){
			col2 = antiAliasing(xp + (cellX * 0.5), yp, cellX * 0.5f, cellY * 0.5f, step);
		}
		if (!(compareVec3(col3, colTot))){
			col3 = antiAliasing(xp, yp + (cellY * 0.5), cellX * 0.5f, cellY * 0.5f, step);
		}
		if (!(compareVec3(col4, colTot))){
			col4 = antiAliasing(xp + (cellX * 0.5), yp + (cellY * 0.5), cellX * 0.5f, cellY * 0.5f, step);
		}
	}
		
	
	return colTot;
}



//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
	float xp, yp;  //grid point
	float cellX = (XMAX-XMIN)/NUMDIV;  //cell width
	float cellY = (YMAX-YMIN)/NUMDIV;  //cell height

	glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glBegin(GL_QUADS);  //Each cell is a quad.

	for(int i = 0; i < NUMDIV; i++)  	//For each grid point xp, yp
	{
		xp = XMIN + i*cellX;
		for(int j = 0; j < NUMDIV; j++)
		{
			
			yp = YMIN + j*cellY;
			glm::vec3 col = antiAliasing(xp, yp, cellX, cellY, 1);
			
			//~ glm::vec3 dir(xp+0.5*cellX, yp+0.5*cellY, -EDIST);	//direction of the primary ray
			//~ glm::vec3 eye(0., 0., 0.);  //The eye position (source of primary rays) is the origin

		    //~ Ray ray = Ray(eye, dir);		//Create a ray originating from the camera in the direction 'dir'
			//~ ray.normalize();				//Normalize the direction of the ray to a unit vector
		    //~ glm::vec3 col = trace (ray, 1); //Trace the primary ray and get the colour value
		    
			glColor3f(col.r, col.g, col.b);
			glVertex2f(xp, yp);				//Draw each cell with its color value
			glVertex2f(xp+cellX, yp);
			glVertex2f(xp+cellX, yp+cellY);
			glVertex2f(xp, yp+cellY);
        }
    }

    glEnd();
    glFlush();
}

void makeCube(){
	
	int x1 = -10.0;
	int x2 = -7.0;
	int y1 = 4.0;
	int y2 = 7.0;
	int z1 = -50;	
	int z2 = -53;
	
	glm::vec3 ptA = glm::vec3(x1, y1, z1);
	glm::vec3 ptB = glm::vec3(x2, y1, z1);
	glm::vec3 ptC = glm::vec3(x2, y1, z2);
	glm::vec3 ptD = glm::vec3(x1, y1, z2);
	glm::vec3 colour = glm::vec3(0.5, 0, 1);
	
	glm::vec3 ptH = glm::vec3(x1, y2, z1);
	glm::vec3 ptG = glm::vec3(x2, y2, z1);
	glm::vec3 ptF = glm::vec3(x2, y2, z2);
	glm::vec3 ptE = glm::vec3(x1, y2, z2);

	Plane *planeBottom = new Plane(ptA, ptB, ptC, ptD, colour);				
	Plane *planeTop= new Plane(ptH, ptG, ptF, ptE, colour);
	Plane *planeFront = new Plane(ptA, ptB, ptG, ptH, colour);				
	Plane *planeBack = new Plane(ptD, ptC, ptF, ptE, colour);
	Plane *planeLeft = new Plane(ptA, ptD, ptE, ptH, colour);
	Plane *planeRight = new Plane(ptB, ptC, ptF, ptG, colour);
	
	sceneObjects.push_back(planeBottom);
	sceneObjects.push_back(planeTop);
	sceneObjects.push_back(planeFront);
	sceneObjects.push_back(planeBack);
	sceneObjects.push_back(planeLeft);
	sceneObjects.push_back(planeRight);
	
}

//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize()
{
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    glClearColor(0, 0, 0, 1);
	texture = TextureBMP("happyElephant.bmp");
	happyTexture = TextureBMP("lemons_pattern-t2.bmp");
	wallTexture = TextureBMP("lemons_pattern-t2.bmp");
	
	//-- Create a pointer to a sphere object
	Sphere *sphere1 = new Sphere(glm::vec3(3.0, 2.0, -125.0), 15.0, glm::vec3(0, 0, 1));
	Sphere *sphere2 = new Sphere(glm::vec3(7.0, -5.0, -80.0), 6.0, glm::vec3(0, 0, 0));
	Sphere *sphere3 = new Sphere(glm::vec3(10.0, 10.0, -60.0), 4.0, glm::vec3(0, 1, 1));
	Sphere *sphere4 = new Sphere(glm::vec3(1-5.0, 15.0, -100.0), 7.0, glm::vec3(0, 1, 1));
	Sphere *sphere5 = new Sphere(glm::vec3(0.0, -12.0, -75.0), 3.0, glm::vec3(0.25, 0.05, 0.5));
	
	Plane *plane = new Plane(glm::vec3(-20., -20, -40),//Point A
							glm::vec3(20., -20, -40),//Point B
							glm::vec3(20., -20, -200),//Point C
							glm::vec3(-20., -20, -200),//Point D
							glm::vec3(1));//Colour
							
	//~ Plane *planeVert = new Plane(glm::vec3(-100., -100, -500),//Point A
						//~ glm::vec3(100., -100, -500),//Point B
						//~ glm::vec3(100., 100, -500),//Point C
						//~ glm::vec3(-100., 100, -500),//Point D
						//~ glm::vec3(1));//Colour
							
	Cylinder *cylinder = new Cylinder(glm::vec3(-15.0, -15.0, -90.0), 6.0, 6.0, glm::vec3(1.0, 1.0, 0));
	
	//--Add the above to the list of scene objects.
	sceneObjects.push_back(sphere1);
	sceneObjects.push_back(sphere2);
	sceneObjects.push_back(sphere3);
	sceneObjects.push_back(plane);
	sceneObjects.push_back(cylinder);
	sceneObjects.push_back(sphere4);
	sceneObjects.push_back(sphere5);

	//sceneObjects.push_back(planeVert);
	
	makeCube();
	
	
}





int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracer");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
