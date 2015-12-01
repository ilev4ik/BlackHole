#pragma once

#include "glm/glm.hpp"
#include "Types.h"
#include "Scene.h"

#include "string"
#include "atlimage.h"

class CTracer
{
public:
	CImage* pImage;
	SRay MakeRay(glm::uvec2 pixelPos, double a, double b);		// cx,cy
	glm::dvec3 TraceRay(SRay, bool, unsigned char*, int, int, int, unsigned char*, int, int, int);			// Trace ray, compute its color
	void RenderImage(const char*, const char*, bool, bool, bool);
	void SaveImageToFile(std::string fileName);
	CImage* LoadImageFromFile(std::string fileName);

	CTracer(double, int, int, double, double, double);
	SCamera m_camera;
	CScene* m_pScene;
};