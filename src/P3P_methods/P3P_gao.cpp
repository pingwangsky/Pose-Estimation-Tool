/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/


#include <gv/P3P_methods/P3P_gao.hpp>
#include <gv/math/roots.hpp>
#include <gv/math/arun.hpp>
#include "iostream"
using namespace std;

// Reference:
/*Gao X S, Hou X R, Tang J, et al. Complete solution classification for
the perspective-three-point problem. Patt Anal Mach Intell,
IEEE Trans[J].IEEE Transactions on Pattern Analysis & Machine Intelligence,
2003, 25(8):930-943.*/
gv::transformations_t
gv::P3P_methods::p3p_gao(const bearingVectors_t & f, const points_t & points)
{
	point_t A = points[0];
	point_t B = points[1];
	point_t C = points[2];

	Eigen::Vector3d tempp;
	tempp = A - B;
	double AB = tempp.norm();
	tempp = B - C;
	double BC = tempp.norm();
	tempp = A - C;
	double AC = tempp.norm();

	bearingVector_t f1 = f[0];
	bearingVector_t f2 = f[1];
	bearingVector_t f3 = f[2];

	double cosalpha = f2.transpose()*f3;
	double cosbeta = f1.transpose()*f3;
	double cosgamma = f1.transpose()*f2;

	double a = pow((BC / AB), 2);
	double b = pow((AC / AB), 2);
	double p = 2 * cosalpha;
	double q = 2 * cosbeta;
	double r = 2 * cosgamma;

	double aSq = a * a;
	double bSq = b * b;
	double pSq = p*p;
	double qSq = q*q;
	double rSq = r*r;

	//if ((pSq + qSq + rSq - p*q*r - 1) == 0)
		//return;

	Eigen::Matrix<double, 5, 1> factors;

	factors[0] = -2 * b + bSq + aSq + 1 - b*rSq*a + 2 * b*a - 2 * a;

	//if (factors[0] == 0)
		//return;

	factors[1] =
		-2 * b*q*a - 2 * aSq*q + b*rSq*q*a - 2 * q + 2 * b*q +
		4 * a*q + p*b*r + b*r*p*a - bSq*r*p;
	factors[2] =
		qSq + bSq*rSq - b*pSq - q*p*b*r + bSq*pSq - b*rSq*a +
		2 - 2 * bSq - a*b*r*p*q + 2 * aSq - 4 * a - 2 * qSq*a + qSq*aSq;
	factors[3] =
		-bSq*r*p + b*r*p*a - 2 * aSq*q + q*pSq*b +
		2 * b*q*a + 4 * a*q + p*b*r - 2 * b*q - 2 * q;
	factors[4] = 1 - 2 * a + 2 * b + bSq - b*pSq + aSq - 2 * b*a;

	std::vector<double> x_temp = math::o4_roots(factors);
	Eigen::Matrix<double, 4, 1> x;
	for (size_t i = 0; i < 4; i++) x[i] = x_temp[i];

	double temp = (pSq*(a - 1 + b) + p*q*r - q*a*r*p + (a - 1 - b)*rSq);
	double b0 = b * temp * temp;

	double rCb = rSq*r;

	Eigen::Matrix<double, 4, 1> tempXP2;
	tempXP2[0] = x[0] * x[0];
	tempXP2[1] = x[1] * x[1];
	tempXP2[2] = x[2] * x[2];
	tempXP2[3] = x[3] * x[3];
	Eigen::Matrix<double, 4, 1> tempXP3;
	tempXP3[0] = tempXP2[0] * x[0];
	tempXP3[1] = tempXP2[1] * x[1];
	tempXP3[2] = tempXP2[2] * x[2];
	tempXP3[3] = tempXP2[3] * x[3];

	Eigen::Matrix<double, 4, 1> ones;
	for (size_t i = 0; i < 4; i++) ones[i] = 1.0;

	Eigen::Matrix<double, 4, 1> b1_part1 =
		(1 - a - b)*tempXP2 + (q*a - q)*x + (1 - a + b)*ones;

	Eigen::Matrix<double, 4, 1> b1_part2 =
		(aSq*rCb + 2 * b*rCb*a - b*rSq*rCb*a - 2 * a*rCb + rCb + bSq*rCb
		- 2 * rCb*b)*tempXP3
		+ (p*rSq + p*aSq*rSq - 2 * b*rCb*q*a + 2 * rCb*b*q - 2 * rCb*q - 2 * p*(a + b)*rSq
		+ rSq*rSq*p*b + 4 * a*rCb*q + b*q*a*rCb*rSq - 2 * rCb*aSq*q + 2 * rSq*p*b*a
		+ bSq*rSq*p - rSq*rSq*p*bSq)*tempXP2
		+ (rCb*qSq + rSq*rCb*bSq + r*pSq*bSq - 4 * a*rCb - 2 * a*rCb*qSq + rCb*qSq*aSq
		+ 2 * aSq*rCb - 2 * bSq*rCb - 2 * pSq*b*r + 4 * p*a*rSq*q + 2 * a*pSq*r*b
		- 2 * a*rSq*q*b*p - 2 * pSq*a*r + r*pSq - b*rSq*rCb*a + 2 * p*rSq*b*q
		+ r*pSq*aSq - 2 * p*q*rSq + 2 * rCb - 2 * rSq*p*aSq*q - rSq*rSq*q*b*p)*x
		+ (4 * a*rCb*q + p*rSq*qSq + 2 * pSq*p*b*a - 4 * p*a*rSq - 2 * rCb*b*q - 2 * pSq*q*r
		- 2 * bSq*rSq*p + rSq*rSq*p*b + 2 * p*aSq*rSq - 2 * rCb*aSq*q - 2 * pSq*p*a
		+ pSq*p*aSq + 2 * p*rSq + pSq*p + 2 * b*rCb*q*a + 2 * q*pSq*b*r + 4 * q*a*r*pSq
		- 2 * p*a*rSq*qSq - 2 * pSq*aSq*r*q + p*aSq*rSq*qSq - 2 * rCb*q - 2 * pSq*p*b
		+ pSq*p*bSq - 2 * pSq*b*r*q*a)*ones;

	Eigen::Matrix<double, 4, 1> b1;
	b1[0] = b1_part1[0] * b1_part2[0];
	b1[1] = b1_part1[1] * b1_part2[1];
	b1[2] = b1_part1[2] * b1_part2[2];
	b1[3] = b1_part1[3] * b1_part2[3];

	Eigen::Matrix<double, 4, 1> y = b1 / b0;
	Eigen::Matrix<double, 4, 1> tempYP2;
	tempYP2[0] = pow(y[0], 2);
	tempYP2[1] = pow(y[1], 2);
	tempYP2[2] = pow(y[2], 2);
	tempYP2[3] = pow(y[3], 2);

	Eigen::Matrix<double, 4, 1> tempXY;
	tempXY[0] = x[0] * y[0];
	tempXY[1] = x[1] * y[1];
	tempXY[2] = x[2] * y[2];
	tempXY[3] = x[3] * y[3];

	Eigen::Matrix<double, 4, 1> v = tempXP2 + tempYP2 - r*tempXY;

	Eigen::Matrix<double, 4, 1> Z;
	Z[0] = AB / sqrt(v[0]);
	Z[1] = AB / sqrt(v[1]);
	Z[2] = AB / sqrt(v[2]);
	Z[3] = AB / sqrt(v[3]);

	Eigen::Matrix<double, 4, 1> X;
	X[0] = x[0] * Z[0];
	X[1] = x[1] * Z[1];
	X[2] = x[2] * Z[2];
	X[3] = x[3] * Z[3];

	Eigen::Matrix<double, 4, 1> Y;
	Y[0] = y[0] * Z[0];
	Y[1] = y[1] * Z[1];
	Y[2] = y[2] * Z[2];
	Y[3] = y[3] * Z[3];

	transformations_t solutions;

	for (int i = 0; i < 4; i++)
	{
		//apply arun to find the transformation
		points_t p_cam;
		p_cam.push_back(X[i] * f1);
		p_cam.push_back(Y[i] * f2);
		p_cam.push_back(Z[i] * f3);

		transformation_t solution = math::arun_complete(points, p_cam);
		//***********************************************
		transformation_t solution1;
		solution1.block<3, 3>(0, 0) = solution.block<3, 3>(0, 0).transpose();
		solution1.col(3) = -solution.block<3, 3>(0, 0).transpose()*solution.col(3);
		//************************************************
		solutions.push_back(solution1);
	}

	return solutions;
}