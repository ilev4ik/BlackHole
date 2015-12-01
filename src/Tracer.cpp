#include <stdlib.h>
#include <stdio.h>
#include "Tracer.h"
#include "glm/gtx/perpendicular.hpp"
#include "glm/gtx/intersect.hpp"
#include <omp.h> 


//atlimage.h for CImage

using namespace glm;

CTracer::CTracer(double alpha, int x, int y, double pos_x, double pos_y, double pos_z)
{
	double rad_alpha = alpha*M_PI / 180.0;
	double d = y / (2 * std::tan(rad_alpha / 2));
	double rad_beta = std::atan2(2 * d, x);

	this->m_camera.m_resolution = glm::uvec2(x, y);
	this->m_camera.m_viewAngle = glm::vec2(rad_alpha, rad_beta);

	this->m_camera.m_pos = glm::dvec3(pos_x, pos_y, pos_z);
	this->m_camera.dist_to_hole = length(this->m_camera.m_pos);

	this->m_camera.m_forward = d * glm::normalize((this->m_camera.m_pos)*(-1.0));
	if (pos_x < 0.0001 && pos_y < 0.0001)
	{
		this->m_camera.m_right = - (double)x * glm::normalize(glm::dvec3(0, 1, 0));
		this->m_camera.m_up = (double)y * glm::normalize(glm::dvec3(1, 0, 0));
	}
	else
	{
		this->m_camera.m_right = (double)x * glm::normalize(glm::dvec3(m_camera.m_forward.y, -m_camera.m_forward.x, 0));
		this->m_camera.m_up = (double)y * glm::normalize(glm::cross(m_camera.m_right, m_camera.m_forward));
	}
}

SRay CTracer::MakeRay(glm::uvec2 pixelPos, double a, double b) //провести луч по (cx,cy) -- экрану проекции, вернуть вектор направления
{
	glm::dvec3 start = m_camera.m_pos;
	glm::dvec3 direction = m_camera.m_forward + (double)((pixelPos.x + a) / m_camera.m_resolution.x - 0.5)*m_camera.m_right +
												(double)((pixelPos.y + b) / m_camera.m_resolution.y - 0.5)*m_camera.m_up;
	return SRay(start, normalize(direction));
}

const double EPS = 0.001;

bool intersect_ray_disk(glm::dvec3 start, glm::dvec3 finish, glm::dvec3 Center, double R_in, double R_out, glm::dvec3& intersect)
{
	double t;
	
	if ((start.z * finish.z <= 0))
	{
		t = -start.z / (finish.z - start.z);
		intersect = dvec3(start.x + t*(finish.x - start.x),
						  start.y + t*(finish.y - start.y),
						  start.z + t*(finish.z - start.z));
		return ((t >= 0) && (t <= 1) && (distance(Center, intersect) <= R_out) && (distance(Center, intersect) >= R_in));
	}
	return false;
}

bool intersect_ray_sphere(glm::dvec3 finish, glm::dvec3 Center, double R)
{
	return (distance(Center, finish) <= R);
}

glm::dvec3 CTracer::TraceRay(SRay ray, bool dop_alpha, 
						     unsigned char* pCurrentLine_disk, int pitch_disk, int H_disk, int W_disk, 
							 unsigned char* pCurrentLine_sky, int pitch_sky, int H_sky, int W_sky)	
{	
	if (pitch_sky < 0)
	{
		pCurrentLine_sky += pitch_sky * (H_sky - 1);
		pitch_sky = -pitch_sky;
	}
	if (pitch_disk < 0)
	{
		pCurrentLine_disk += pitch_disk * (H_disk - 1);
		pitch_disk = -pitch_disk;
	}
	unsigned char r, g, b, alpha = 0;
	glm::dvec3 color(1, 1, 1), color_0;
	glm::dvec3 pos1 = ray.m_start, pos0;				
	glm::dvec3 V = lightSpeed * ray.m_dir;
	double GM = Gravy * this->m_pScene->getM();
	double dt = 5;
	glm::dvec3 a, an;
	double R_hole = this->m_pScene->getRH();
	double R_disk = this->m_pScene->getRD();
	glm::dvec3 Center = this->m_pScene->getCenter();

	bool found = false;
	glm::dvec3 intersect_disk_point;
	int iterations = 0;
	bool sphere = false, disk = false;

	int max_iter;
	if (ray.m_start.x == 0 && ray.m_start.y == 0)
		max_iter = round(abs(2 * ray.m_start.z / (dt * lightSpeed)));
	else
		max_iter = round((4 * R_hole + length(dvec3(ray.m_start.x, ray.m_start.y, 0))) / (dt * lightSpeed));

	while (!found && iterations < max_iter)
	{
		iterations++;
		a = -GM * pos1 / pow(length(pos1), 3);
		an = perp(a, V);
		V += an*dt;
		V = lightSpeed*normalize(V);
		pos0 = pos1;
		pos1 += dt*V + dt*dt*0.5*an;
		sphere = intersect_ray_sphere(pos1, Center, R_hole);
		disk = intersect_ray_disk(pos0, pos1, Center, R_hole, R_disk, intersect_disk_point);

		if (sphere)
		{
			color = dvec3(0, 0, 0);
			found = true;
		}
		else
			if (disk)
			{
				int j = round(W_disk * 0.5 * (intersect_disk_point.x / R_disk + 1));
				int i = round(H_disk * 0.5 * (intersect_disk_point.y / R_disk + 1));
				if (i >= H_disk)
					i = H_disk - 1;
				if (j >= W_disk)
					j = W_disk - 1;

				b = pCurrentLine_disk[i*pitch_disk + j * 4];
				g = pCurrentLine_disk[i*pitch_disk + j * 4 + 1];
				r = pCurrentLine_disk[i*pitch_disk + j * 4 + 2];
				alpha = pCurrentLine_disk[i*pitch_disk + j * 4 + 3];

				if (alpha != 0)
				{
					color_0 = dvec3(r, g, b);
					if (!dop_alpha)
					{
						found = true;
						color = color_0;
					}
				}
			}
	}
	//значит фон
	if (iterations == max_iter)
	{
		pos1 = normalize(pos1 - pos0);
		double phi = atan2(pos1.x, pos1.y);
		double tetta = asin(pos1.z);

		int j = round(W_sky / 2 * (phi / M_PI + 1));
		int i = round(H_sky / 2 * (2 * tetta / M_PI + 1));
		
		b = pCurrentLine_sky[i*pitch_sky + j * 3];
		g = pCurrentLine_sky[i*pitch_sky + j * 3 + 1];
		r = pCurrentLine_sky[i*pitch_sky + j * 3 + 2];
		color = dvec3(r, g, b);
	}

	if (dop_alpha)
		color = (double)alpha / 255 * color_0 + (double)(1 - alpha / 255)*color;

	return (1/255.0)*color;
}
	

void CTracer::RenderImage(const char* background, const char* disk, bool dop_alpha, bool dop_openmp, bool dop_aliasing)
{
		// Rendering
	size_t xRes = m_camera.m_resolution.x;
	size_t yRes = m_camera.m_resolution.y;

	CImage* pImage_S = LoadImageFromFile(background);
	auto pData_S = (unsigned char*)pImage_S->GetBits();
	auto pCurrentLine_sky = pData_S;
	int pitch_S = pImage_S->GetPitch();
	int H_sky = pImage_S->GetHeight();
	int W_sky = pImage_S->GetWidth();
	
	CImage* pImage_D = LoadImageFromFile(disk);
	auto pData_D = (unsigned char*)pImage_D->GetBits();
	auto pCurrentLine_D = pData_D;
	int pitch_D = pImage_D->GetPitch();
	int H_disk = pImage_D->GetHeight();
	int W_disk = pImage_D->GetWidth();

	this->m_camera.m_pixels.resize(xRes*yRes);

	if (dop_openmp)
	{
		for (int i = 0; i < yRes; i++)
		{
			#pragma omp parallel for
			for (int j = 0; j < xRes; j++)
			{
				SRay ray = MakeRay(uvec2(j, i), 0.5, 0.5);
				glm::dvec3 result = TraceRay(ray, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky);

				if (dop_aliasing)
				{
					SRay ray00 = MakeRay(uvec2(j, i), 0, 0);
					SRay ray01 = MakeRay(uvec2(j, i), 0, 1);
					SRay ray10 = MakeRay(uvec2(j, i), 1, 0);
					SRay ray11 = MakeRay(uvec2(j, i), 1, 1);

					glm::dvec3 color = TraceRay(ray00, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky)
								     + TraceRay(ray01, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky)
								     + TraceRay(ray10, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky)
						             + TraceRay(ray11, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky);
					result = (1 / 5.0)*(color + result);
				}
				m_camera.m_pixels[i * xRes + j] = result;
				
			}
			printf("(%d)\n", i);
		}
	}
	else
		for (int i = 0; i < yRes; i++)
		{
			for (int j = 0; j < xRes; j++)
			{
				SRay ray = MakeRay(uvec2(j, i), 0.5, 0.5);
				glm::dvec3 result = TraceRay(ray, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky);

				if (dop_aliasing)
				{
					SRay ray00 = MakeRay(uvec2(j, i), 0, 0);
					SRay ray01 = MakeRay(uvec2(j, i), 0, 1);
					SRay ray10 = MakeRay(uvec2(j, i), 1, 0);
					SRay ray11 = MakeRay(uvec2(j, i), 1, 1);

					glm::dvec3 color = TraceRay(ray00, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky)
						+ TraceRay(ray01, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky)
						+ TraceRay(ray10, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky)
						+ TraceRay(ray11, dop_alpha, pCurrentLine_D, pitch_D, H_disk, W_disk, pCurrentLine_sky, pitch_S, H_sky, W_sky);
					result = (1 / 5.0)*(color + result);
				}
				
				m_camera.m_pixels[i * xRes + j] = result;
			}
			printf("(%d)\n", i);
		}
}

void CTracer::SaveImageToFile(std::string fileName)
{
  CImage image;

  int width = m_camera.m_resolution.x;
  int height = m_camera.m_resolution.y;

  image.Create(width, height, 24);
    
	int pitch = image.GetPitch();
	unsigned char* imageBuffer = (unsigned char*)image.GetBits();

	if (pitch < 0)
	{
		imageBuffer += pitch * (height - 1);
		pitch =- pitch;
	}

	int i, j;
	int imageDisplacement = 0;
	int textureDisplacement = 0;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
		  vec3 color = m_camera.m_pixels[textureDisplacement + j];

		  imageBuffer[imageDisplacement + j * 3] = clamp(color.b, 0.0f, 1.0f) * 255.0f;
		  imageBuffer[imageDisplacement + j * 3 + 1] = clamp(color.g, 0.0f, 1.0f) * 255.0f;
		  imageBuffer[imageDisplacement + j * 3 + 2] = clamp(color.r, 0.0f, 1.0f) * 255.0f;
		}

		imageDisplacement += pitch;
		textureDisplacement += width;
	}

	image.Save(fileName.c_str());
	image.Destroy();
}

CImage* CTracer::LoadImageFromFile(std::string fileName)
{
  CImage* pImage = new CImage;

  if(SUCCEEDED(pImage->Load(fileName.c_str())))
    return pImage;
  else
  {
    delete pImage;
    return NULL;
  }
}