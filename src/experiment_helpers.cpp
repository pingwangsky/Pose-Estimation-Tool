/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/


#include "gv\experiment_helpers.hpp"

//*******************MY WORK BEGIN**********************
//my work
//自己添加的程序
void
gv::GeneratePnPdata(
points_t & worldpoints,
Image_points & image_points,
translation_t & position,
rotation_t & rotation,
//bearingVectors_t & bearingVectors,
size_t numberPoints,
double depth,
double noise)
{
	double focal = 800;
	points_t camerapoints;
	point_t point;

	srand((unsigned)time(NULL));

	for (size_t i = 0; i<numberPoints; i++)
	{
		point[0] = -4 + (double)rand() / ((double)RAND_MAX)*(4 + 4);
		point[1] = -4 + (double)rand() / ((double)RAND_MAX)*(4 + 4);
		point[2] = depth + (double)rand() / ((double)RAND_MAX) * 8;

		camerapoints.push_back(point);
	}

	point_t mean_t(0, 0, 0);
	for (size_t i = 0; i < numberPoints; i++)
	{
		mean_t[0] = (double)mean_t[0] + camerapoints[i][0];
		mean_t[1] = (double)mean_t[1] + camerapoints[i][1];
		mean_t[2] = (double)mean_t[2] + camerapoints[i][2];
	}
	mean_t = mean_t / numberPoints;
	position = mean_t;		//位置矩阵
	rotation = generateRandomRotation(0.5);		//旋转矩阵 
	point_t point1;
	for (size_t i = 0; i<numberPoints; i++)
	{
		point1 = rotation.transpose()*(camerapoints[i] - mean_t);
		worldpoints.push_back(point1);
	}
	Image_point image_point;

	for (size_t i = 0; i<numberPoints; i++)
	{
		image_point[0] = (camerapoints[i][0] / camerapoints[i][2] * focal + (double)rand() / ((double)RAND_MAX)*noise) / focal;
		image_point[1] = (camerapoints[i][1] / camerapoints[i][2] * focal + (double)rand() / ((double)RAND_MAX)*noise) / focal;
		image_points.push_back(image_point);
	}

	point_t temp;
	for (size_t i = 0; i<numberPoints; i++)
	{
		temp[0] = image_points[i][0];
		temp[1] = image_points[i][1];
		temp[2] = 1;
		//bearingImage.push_back(temp);
		//temp = temp / temp.norm();
		//bearingVectors.push_back(temp);
	}
}

//show experiment data;

//用来显示实验数据
void
gv::ShowExperimentData(
points_t & worldpoints,
Image_points & image_points,
const translation_t & position,
const rotation_t & rotation,
/*bearingVectors_t & bearingVectors,*/
size_t numberPoints,
double depth,
double noise)
{
	Eigen::MatrixXd show_worldpoints(3,numberPoints);
	Eigen::MatrixXd show_image_points(2, numberPoints);
	Eigen::MatrixXd show_bearingVectors(3,numberPoints);

	for (size_t i = 0; i<numberPoints; i++)
	{
		show_worldpoints.col(i) = worldpoints[i];

	}

	for (size_t i = 0; i<numberPoints; i++)
	{
		show_image_points.col(i) = image_points[i];
	}


	/*for (size_t i = 0; i<numberPoints; i++)
	{
		show_bearingVectors.col(i) = bearingVectors[i];
	}*/



	std::cout << "world points is:" << std::endl;
	std::cout << show_worldpoints << std::endl << std::endl;
	std::cout << "image points is:" << std::endl;
	std::cout << show_image_points << std::endl << std::endl;
	//std::cout << "bearingVectors is:" << std::endl;
	//std::cout << show_bearingVectors << std::endl << std::endl;
	std::cout << "the random position is:" << std::endl;
	std::cout << position << std::endl << std::endl;
	std::cout << "the random rotation is:" << std::endl;
	std::cout << rotation << std::endl << std::endl;
	std::cout << "the noise in the data is:" << std::endl;
	std::cout << noise << std::endl;
	std::cout << "the depth in the data is:" << std::endl;
	std::cout << depth << std::endl;
}

//显示最终解
void
gv::ShowSolutions(transformations_t & solutions)
{
	for (size_t i = 0; i < solutions.size(); i++)
		std::cout << solutions[i] << std::endl << std::endl;
}

//***************我添加的程序结束*********************
//*******************MY WORK END**********************

