/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/

#ifndef GV_SRPnP_MAIN_HPP_
#define GV_SRPnP_MAIN_HPP_

#include <stdlib.h>
#include <gv/types.hpp>

namespace gv
{
namespace SRPnP
{
	//My simple, robust and fast method(SRPnP);
	gv::transformation_t srpnp(const Image_points & xn, const points_t & Xw);

	Eigen::MatrixXd getp3p(double l1, double l2, double A5, double C1, double C2, double D1, double D2, double D3);

	Eigen::MatrixXd getpoly7(Eigen::MatrixXd F);

	Eigen::MatrixXd fucntionB(rotation_t R, point_t P, double ui, double vi);

	Eigen::MatrixXd MatrixCs(point_t P, double ui, double vi);

	std::vector<double> Filtersolution(complex_C c);

	Eigen::Vector4d matrix2quaternion(rotation_t R);

	Eigen::Vector3d RefineGaussNewton(Eigen::Vector3d solution, Eigen::MatrixXd E, Eigen::MatrixXd G1);

	Eigen::Vector3d Cayley(rotation_t R);

}
}

#endif 
