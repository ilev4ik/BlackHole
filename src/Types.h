#pragma once

#define _USE_MATH_DEFINES
#include "glm/glm.hpp" 
#include "vector"
#include <cmath>
#include <ctime>

class SRay
{
public:
	glm::dvec3 m_start; //������
	glm::dvec3 m_dir;	//����
public:
	SRay(glm::dvec3 a, glm::dvec3 b) : m_start(a), m_dir(b) {};

};

struct SCamera
{
	glm::dvec3 m_pos;          // ������� ������
	
	glm::dvec3 m_forward;      
	glm::dvec3 m_right;
	glm::dvec3 m_up;
	double dist_to_hole;
	glm::vec2 m_viewAngle;    // ���� ������ � �������� (�,�)
	glm::uvec2 m_resolution;  // ����������
	std::vector<glm::vec3> m_pixels;  // ������ �������� �������� �� �����
};

//�����
struct SMesh
{
	std::vector<glm::dvec4> m_vertices;  // ���������� ������
};
