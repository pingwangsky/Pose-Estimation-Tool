/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/


#include <gv/P3P_methods/P3P_li.hpp>
#include <gv/math/roots.hpp>
#include <gv/math/arun.hpp>
#include "iostream"
using namespace std;

// Author: Ping Wang
// Li's PST p3p solver;
// Reference:
/* SHIQILI, CHIXU. A STABLE DIRECT SOLUTION OF
* PERSPECTIVE-THREE-POINT PROBLEM[J]. International
* Journal of Pattern Recognition & Artificial Intelligence,
* 2011, 25(5):627-642.*/
gv::transformations_t
gv::P3P_methods::p3p_li(const bearingVectors_t & f, const points_t & p)
{
	point_t P1 = p[0];
	point_t P2 = p[1];
	point_t P3 = p[2];

	bearingVector_t f1 = f[0];
	bearingVector_t f2 = f[1];
	bearingVector_t f3 = f[2];

	double D1 = (P2 - P1).norm();
	double D2 = (P3 - P1).norm();
	double D3 = (P3 - P2).norm();

	double cg1 = f1.transpose()*f2;
	double cg2 = f1.transpose()*f3;
	double cg3 = f2.transpose()*f3;
	double sg1 = sqrt(1 - pow(cg1, 2));
	double sg2 = sqrt(1 - pow(cg2, 2));
	double sg3 = sqrt(1 - pow(cg3, 2));

	double s0 = 1;
	double s1 = s0*cg1;
	double s2 = s0*cg2;

	double C1 = s0*sg1;
	double C2 = s0*sg2;

	double A1 = pow(D2 / D1, 2);
	double A2 = A1*pow(C1, 2) - pow(C2, 2);
	double A3 = s2*cg3 - s1;
	double A4 = s1*cg3 - s2;
	double A5 = cg3;
	double A6 = (pow(D3, 2) - pow(D1, 2) - pow(D2, 2)) / (2 * pow(D1, 2));
	double A7 = pow(s0, 2) - pow(s1, 2) - pow(s2, 2) + s1*s2*cg3 + A6*pow(C1, 2);

	Eigen::Matrix<double, 5, 1> factors;

	factors(0, 0) = pow(A6, 2) - A1*pow(A5, 2);
	factors(1, 0) = 2 * (A3*A6 - A1*A4*A5);
	factors(2, 0) = pow(A3, 2) + 2 * A6*A7 - A1*pow(A4, 2) - A2*pow(A5, 2);
	factors(3, 0) = 2 * (A3*A7 - A2*A4*A5);
	factors(4, 0) = pow(A7, 2) - A2*pow(A4, 2);

	std::vector<double> realRoots = math::o4_roots(factors);

	transformations_t solutions;

	//retrieve orientation and position;
	for (int i = 0; i < 4; i++)
	{
		double t1 = realRoots[i];
		double c = A4 + A5*t1;
		double t2 = -(A3*t1 + A6*pow(t1, 2) + A7) / (A4 + A5*t1);
		double lambda = D1 / sqrt(pow(t1, 2) + pow(C1, 2));
		double d1 = lambda*s0;
		double d2 = lambda*(s1 + t1);
		double d3 = lambda*(s2 + t2);

		//apply arun to find the transformation
		points_t p_cam;
		p_cam.push_back(f1*d1);
		p_cam.push_back(f2*d2);
		p_cam.push_back(f3*d3);
		transformation_t solution = math::arun_complete(p, p_cam);
		transformation_t solution1;
		solution1.block<3, 3>(0, 0) = solution.block<3, 3>(0, 0).transpose();
		solution1.col(3) = -solution.block<3, 3>(0, 0).transpose()*solution.col(3);
		solutions.push_back(solution1);
	}

	return solutions;
}
