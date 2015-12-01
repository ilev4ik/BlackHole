#pragma once

#include "Types.h"
const double	lightSpeed = 3e+8, Gravy = 6.674e-11;

class CScene
{
public:
//	SMesh Hole;
//	SMesh Disk;
	glm::dvec3 holeCenter;
	double RadHole, RadDisk;					
	unsigned int xRes, yRes;
	double Mass;
public:
	CScene(double, double);
//	glm::dmat4x4 rot(double);
//	glm::dmat4x4 par(double, double, double);
	glm::dvec3 getCenter();
	double getRH();
	double getRD();
	double getM();
//	void buildHole();
//	std::vector<glm::dvec4> cirlce(glm::vec3, double);
//	void buildDisk();
};