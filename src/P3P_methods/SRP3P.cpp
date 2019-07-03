/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/


#include <gv/P3P_methods/SRP3P.hpp>
#include <gv/math/roots.hpp>
#include "iostream"
using namespace std;


// Author: Ping Wang
// My simple and robust p3p solver;
gv::transformations_t
gv::P3P_methods::srp3p_wang(const bearingVectors_t & f, const points_t & p)
{
	point_t P1 = p[0];
    point_t P2 = p[1];
    point_t P3 = p[2];

	bearingVector_t f1 = f[0];
    bearingVector_t f2 = f[1];
    bearingVector_t f3 = f[2];

	double N1 = pow((P2 - P1).norm(),2);
	double N2 = pow((P3 - P1).norm(),2);
	Eigen::Vector3d L1 = P3 - P1;
	Eigen::Vector3d L2 = P2 - P1;
	double N3=L1(0)*L2(0)+L1(1)*L2(1)+L1(2)*L2(2);

	double cos1=f1.transpose()*f2;
	double cos2=f1.transpose()*f3;
	double cos3=f2.transpose()*f3;


	double sin1=sqrt(1-pow(cos1,2));
	double sin2=sqrt(1-pow(cos2,2));

	double A=pow(sin1,2);
	double B=pow(sin2,2);
	double C1=cos2*cos3-cos1;
	double C2=cos1*cos3-cos2;
	double C3=cos3;
	double C4=1+cos1*cos2*cos3-pow(cos1,2)-pow(cos2,2);

	double K1=N1/N2;
	double K2=N1/N3;

	double D=A-K1*B;
	double E1=K2*C1;
	double E2=K2*C2;
	double E3=K2*C3;
	double E4=A-K2*C4;

	Eigen::Matrix<double,5,1> factors;

	factors(0,0)=pow(E3,2)-K1;
	factors(1,0)=2*E2*E3+2*K1*E1;
	factors(2,0)=pow(E2,2)-2*E4*K1-K1*pow(E1,2)+pow(E3,2)*D;
	factors(3,0)=2*K1*E1*E4+2*E2*E3*D;
	factors(4,0)=pow(E2,2)*D-K1*pow(E4,2);

	std::vector<double> realRoots = math::o4_roots(factors);

    transformations_t solutions;

	//retrieve the orientation and position;
	for( int i = 0; i < 4; i++ )
    {
		double x=realRoots[i];
		double y=(pow(x,2)-x*E1+E4)/(x*E3+E2);
		double lambda1=x+cos1;
		double lambda2=y+cos2;
		double d1=sqrt(N1/(pow(x,2)+A));
		double d2=lambda1*d1;
		double d3=lambda2*d1;
		point_t Xc1=f1*d1;
		point_t Xc2=f2*d2;
		point_t Xc3=f3*d3;

	    //cout<<realRoots[3] <<endl<<endl;
		
		Eigen::Vector3d m1=(P2-P1)/(P2-P1).norm();
		Eigen::Vector3d m2=(P3-P1)/(P3-P1).norm();
		Eigen::Vector3d m3=m1.cross(m2);
		m3=m3/m3.norm();
		rotation_t m;
		m.col(0) = m1;
        m.col(1) = m2;
        m.col(2) = m3;

		Eigen::Vector3d n1=(Xc2-Xc1)/(Xc2-Xc1).norm();
		Eigen::Vector3d n2=(Xc3-Xc1)/(Xc3-Xc1).norm();
		Eigen::Vector3d n3=n1.cross(n2);
		n3=n3/n3.norm();
		rotation_t n;
		n.col(0) = n1;
        n.col(1) = n2;
        n.col(2) = n3;

        transformation_t solution;		//3*4ÁÐµÄ¾ØÕó
        solution.block<3,3>(0,0) = n*m.inverse();
		solution.col(3) = Xc1-n*m.inverse()*P1;

        solutions.push_back(solution);
    }

	return solutions;
}
