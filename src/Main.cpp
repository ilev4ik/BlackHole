#include "Tracer.h"
#include "stdio.h"


int main(int argc, char** argv)
{
//	printf("Time: ");
//	clock_t time = clock();

	int xRes = 1920;
	int yRes = 1200;
	double M = 8.57e+36, alpha = 60, pos_x = 12.0e+10, pos_y = 12.0e+10, pos_z = 1.0e+10;
	int c = 3;
	char background[256], disk[256];
	char c_alpha = 'T', c_openmp = 'T', c_anti_aliasing = 'F';

	if(argc == 2) // There is input file in parameters
	{
		FILE* file = fopen(argv[1], "r");
		if(file)
		{
			if (fscanf(file, "Background texture path: %s\n Disk texture path: %s\n Resolution: %d %d\n Hole mass: %lf\n Ratio: %d\n View angle: %lf\n Camera position: %lf %lf %lf\n",
				background, disk, &xRes, &yRes, &M, &c, &alpha, &pos_x, &pos_y, &pos_z) != 10)
			{
				printf("Invalid config format! Using default parameters.\r\n");
			}
			if (fscanf(file, "Alpha blending: %c\n OpenMP: %c\n Anti-aliasing: %c\n", &c_alpha, &c_openmp, &c_anti_aliasing) != 3)
			{
				printf("Invalid config format! Using default parameters.\r\n");
			}
			fclose(file);
		}
		else
			printf("Invalid config path! Using default parameters.\r\n");
	}
	else
		printf("No config! Using default parameters.\r\n");


	bool dop_alpha = c_alpha == 'T';
	bool dop_openmp = c_openmp == 'T';
	bool dop_aliasing = c_anti_aliasing == 'T';
	CTracer tracer(alpha, xRes, yRes, pos_x, pos_y, pos_z);
	CScene scene(M,c);

	tracer.m_pScene = &scene;
	//tracer.m_pScene->buildDisk();
	//tracer.m_pScene->buildHole();
	tracer.RenderImage(background, disk, dop_alpha, dop_openmp, dop_aliasing);
	tracer.SaveImageToFile("Result.png");

//	printf("%lf mins.\n", (double)(clock() - time) / CLOCKS_PER_SEC / 60);

//	system("Pause");
	return 0;
}