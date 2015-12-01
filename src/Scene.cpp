#include "Scene.h"
extern unsigned int W,H;

glm::dvec4 operator* (glm::dmat4x4 mat, glm::dvec4 vec)
{
	glm::dvec4 result;
	for (int i = 0; i < 4; ++i)
	{
		result[i] = glm::dot(mat[i], vec);
	}

	return result;
}

double F(double a, double b, double R)
{
	return (a*a + b*b - R*R);
}

CScene::CScene(double M, double ratio)
{
	this->Mass = M;
	this->RadHole = std::round(((2.0*Gravy*Mass / lightSpeed) / lightSpeed));
	this->RadDisk = ratio*RadHole;
	holeCenter = glm::dvec3(0.0, 0.0, 0.0);

}

glm::dvec3 CScene::getCenter()
{
	return this->holeCenter;
}

double CScene::getM()
{
	return this->Mass;
}

double CScene::getRH()
{
	return this->RadHole;
}

double CScene::getRD()
{
	return this->RadDisk;
}
/*
std::vector <glm::dvec4> CScene::cirlce(glm::vec3 O, double R)
{
	double Y = O.y;
	std::vector <glm::dvec4> buffer1;
	glm::dvec4 M = glm::dvec4(O.x, O.y + R, O.z, 1);
	buffer1.push_back(M);
	
	do
	{
		M.x += 1;
		double y1 = M.y;
		double y2 = M.y - 1;
		
		while (F(M.x - O.x, y1 - O.y, R) > 0 && F(M.x - O.x, y2 - O.y, R) > 0)
		{
			y1 = y2;
			--y2;
		}
		
		M.y = y1 - 0.5;
		double Fm = F(M.x - O.x, M.y - O.y, R);
		if (Fm < 0)
		{
			M.y += 0.5;
			glm::dvec4 vert = glm::dvec4(M.x, M.y, M.z, 1);
			buffer1.push_back(vert);
		}
		else if (Fm >= 0)
		{
			M.y -= 0.5;
			glm::dvec4 vert = glm::dvec4(M.x, M.y, M.z, 1);
			buffer1.push_back(vert);
		}
	} while (M.y != O.y);

//----------------------buffer2----------------------------------

	std::vector <glm::dvec4> buffer2;
	int size = buffer1.size();
	for (int i = 0; i < size; ++i)
	{
		glm::dvec4 temp = par(-O.x, -O.y, -O.z)*buffer1[i];	//перенос в ноль
		temp = rot(M_PI / 2)*temp;							//поворот относительно нуля на 90
		temp = par(O.x, O.y, O.z)*temp;						//обратный перенос на место :)
		buffer2.push_back(temp);
		
	}
	buffer2.erase(buffer2.begin());							//вынимаем левую точку,
	buffer2.pop_back();										//вынимаем правую точку для корректности 
	
	//concat
	buffer1.insert(buffer1.end(), buffer2.begin(), buffer2.end());
	buffer2.clear();											
	
	size = buffer1.size();									//очищенный буфер для дальнейшего использования
	for (int i = 0; i < size; ++i)
	{	
		glm::dvec4 temp = par(-O.x, -O.y, -O.z)*buffer1[i];	
		temp = rot(M_PI)*temp;						
		temp = par(O.x, O.y, O.z)*temp;					
		buffer1.push_back(temp);
	}

	return buffer1;
}
*/
/*
glm::dmat4x4 CScene::rot(double alpha)
{
	glm::dmat4x4 rot;
	rot[0] = glm::dvec4(cos(alpha), -sin(alpha), 0, 0);
	rot[1] = glm::dvec4(sin(alpha), cos(alpha), 0, 0);
	rot[2] = glm::dvec4(0, 0, 1, 0);
	rot[3] = glm::dvec4(0, 0, 0, 1);
	return rot;
}

glm::dmat4x4 CScene::par(double x, double y, double z)
{
	glm::dmat4x4 parxy;
	parxy[0] = glm::dvec4(1, 0, 0, x);
	parxy[1] = glm::dvec4(0, 1, 0, y);
	parxy[2] = glm::dvec4(0, 0, 1, z);
	parxy[3] = glm::dvec4(0, 0, 0, 1);
	return parxy;
}

void CScene::buildDisk()
{
	clock_t time = clock();
	this->Disk.m_vertices = cirlce(holeCenter, RadHole);
	printf("Disk building time: %lf\n", double(clock() - time) / CLOCKS_PER_SEC);
	return;
}

void CScene::buildHole()
{
	clock_t time = clock();
	glm::vec3 centerPos_up = this->holeCenter;
	glm::vec3 centerPos_down = this->holeCenter;

	for (int z_Step = 0; z_Step < this->RadHole; ++z_Step)
	{
		std::vector <glm::dvec4> temp_up = cirlce(centerPos_up, this->RadHole - z_Step);	//z_Step ~ RadHole
		this->Hole.m_vertices.insert(this->Hole.m_vertices.end(), temp_up.begin(), temp_up.end());

		std::vector <glm::dvec4> temp_down = cirlce(centerPos_down, this->RadHole - z_Step);	//z_Step ~ RadHole
		this->Hole.m_vertices.insert(this->Hole.m_vertices.end(), temp_down.begin(), temp_down.end());

		(centerPos_up.z)++;
		(centerPos_down.z)--;
	}
	printf("Hole building time: %lf\n", double(clock() - time) / CLOCKS_PER_SEC);
}
*/
