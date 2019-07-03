/******************************************************************************
* This programe is implemented in Visual Studio 2015.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/


#include <gv/P3P_methods/EP3P.hpp>
#include <gv/math/roots.hpp>
#include "iostream"
using namespace std;

// Author: Ping Wang
// My efficient p3p solver;
gv::transformations_t
gv::P3P_methods::ep3p_wang(const bearingVectors_t & f, const points_t & p)
{
	point_t P1 = p[0];
	point_t P2 = p[1];
	point_t P3 = p[2];

	bearingVector_t f1 = f[0];
	bearingVector_t f2 = f[1];
	bearingVector_t f3 = f[2];

	Eigen::Vector3d uz = (P2 - P1) / (P2 - P1).norm();
	Eigen::Vector3d uy = uz.cross(P3 - P1) / uz.cross(P3 - P1).norm();
	Eigen::Vector3d ux = uy.cross(uz);

	rotation_t R_temp;
	R_temp.row(0) = ux.transpose();
	R_temp.row(1) = uy.transpose();
	R_temp.row(2) = uz.transpose();

	point_t P_u1(0, 0, 0);
	point_t P_u2 = R_temp*(P2 - P1);
	point_t P_u3 = R_temp*(P3 - P1);

	double z_2 = P_u2[2];
	double x_3 = P_u3[0];
	double z_3 = P_u3[2];

	//-----------------------
	double A1 = f1[0] / f1[2];
	double A2 = f1[1] / f1[2];
	double B1 = f2[0] / f2[2];
	double B2 = f2[1] / f2[2];
	double C1 = f3[0] / f3[2];
	double C2 = f3[1] / f3[2];
	//-----------------------
	double T1 = (B1 - A1) / z_2;
	double T2 = (B2 - A2) / z_2;
	//-----------------------
	double D1 = (C1*z_3 - B1*z_3) / x_3;
	double D2 = (C1 - A1 - z_3*T1) / x_3;
	double D3 = (C2*z_3 - B2*z_3) / x_3;
	double D4 = (C2 - A2 - z_3*T2) / x_3;
	//-----------------------
	double G1 = D1*B1 + D3*B2;
	double G2 = C1*B1 + C2*B2 + 1;
	double G3 = B1*D2 + D1*T1 + D3*T2 + B2*D4;
	double G4 = C1*T1 + C2*T2;
	double G5 = D2*T1 + D4*T2;
	//-----------------------
	double H1 = pow(C1, 2) + pow(C2, 2) + 1;
	double H2 = pow(D1, 2) + pow(D3, 2) - pow(B1, 2) - pow(B2, 2) - 1;
	double H3 = 2 * C1*D1 + 2 * C2*D3;
	double H4 = 2 * C1*D2 + 2 * C2*D4;
	double H5 = 2 * D1*D2 + 2 * D3*D4 - 2 * B1*T1 - 2 * B2*T2;
	double H6 = pow(D2, 2) + pow(D4, 2) - pow(T1, 2) - pow(T2, 2);
	//-----------------------
	Eigen::Matrix<double, 5, 1> factors;
	factors(0, 0) = H1*pow(G1, 2) + H2*pow(G2, 2) - H3*G1*G2;
	factors(1, 0) = 2 * H1*G1*G3 + 2 * H2*G2*G4 - H3*G1*G4 - H3*G2*G3
		- H4*G1*G2 + H5*pow(G2, 2);
	factors(2, 0) = H1*pow(G3, 2) + 2 * H1*G1*G5 + H2*pow(G4, 2) - H3*G3*G4
		- H3*G2*G5 - H4*G1*G4 - H4*G2*G3 + 2 * H5*G2*G4 + H6*pow(G2, 2);
	factors(3, 0) = 2 * H1*G3*G5 - H3*G4*G5 - H4*G3*G4 - H4*G2*G5 +
		H5*pow(G4, 2) + 2 * H6*G2*G4;
	factors(4, 0) = H1*pow(G5, 2) - H4*G4*G5 + H6*pow(G4, 2);

	std::vector<double> realRoots = math::o4_roots(factors);

	transformations_t solutions;

	for (int i = 0; i < realRoots.size(); i++)
	{
		double S4 = realRoots[i];
		double S7 = -(G1*pow(S4, 2) + G3*S4 + G5) / (G2*S4 + G4);
		double S1 = A1;
		double S2 = A2;
		double S3 = B1*S4 + T1;
		double S5 = B2*S4 + T2;
		double S6 = C1*S7 + D1*S4 + D2;
		double S8 = C2*S7 + D3*S4 + D4;
		double tz = 1 / sqrt(pow(S3, 2) + pow(S4, 2) + pow(S5, 2));
		double tx = S1*tz;
		double ty = S2*tz;
		double r1 = S6*tz;
		double r3 = S3*tz;
		double r4 = S8*tz;
		double r6 = S5*tz;
		double r7 = S7*tz;
		double r9 = S4*tz;
		translation_t T_F(tx, ty, tz);
		double r2 = r6*r7 - r4*r9;
		double r5 = r1*r9 - r3*r7;
		double r8 = r3*r4 - r1*r6;
		rotation_t R_F;
		R_F << r1, r2, r3, r4, r5, r6, r7, r8, r9;

		transformation_t solution;		//3*4 matrix;
		solution.block<3, 3>(0, 0) = R_F*R_temp;
		solution.col(3) = T_F - R_F*R_temp*P1;

		solutions.push_back(solution);
	}

	return solutions;

}

/*std::vector<double>
gv::SRPnP::Filtersolution(complex_C c)
{
//select minima;
double maxreal = abs(c[0].real());
for (unsigned int i = 0; i < c.size(); i++)
{
if (abs(c[i].real())>maxreal)
{
maxreal = abs(c[i].real());
}
}

std::vector<double> realRoots;
for (unsigned int i = 0; i < c.size(); i++)
{
if (abs(c[i].imag()) / maxreal<0.001)
{
realRoots.push_back(c[i].real());
}
}

return realRoots;
}*/