#pragma once

#define _USE_MATH_DEFINES
#include "glm/glm.hpp" 
#include "vector"
#include <cmath>
#include <ctime>

class SRay
{
public:
	glm::dvec3 m_start; //откуда
	glm::dvec3 m_dir;	//куда
public:
	SRay(glm::dvec3 a, glm::dvec3 b) : m_start(a), m_dir(b) {};

};

struct SCamera
{
	glm::dvec3 m_pos;          // Позиция камеры
	
	glm::dvec3 m_forward;      
	glm::dvec3 m_right;
	glm::dvec3 m_up;
	double dist_to_hole;
	glm::vec2 m_viewAngle;    // Углы обзора в радианах (х,у)
	glm::uvec2 m_resolution;  // Разрешение
	std::vector<glm::vec3> m_pixels;  // Массив пикселей проекции на линзу
};

//Сетка
struct SMesh
{
	std::vector<glm::dvec4> m_vertices;  // координаты вершин
};
