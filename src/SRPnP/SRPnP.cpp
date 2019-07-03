/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/



#include <gv/SRPnP/SRPnP.hpp>
#include <gv/math/roots.hpp>
#include <gv/math/arun.hpp>
#include "iostream"
using namespace std;

//my simple, robust and fast method for PnP problem;
gv::transformation_t
gv::SRPnP::srpnp(const Image_points & xn, const points_t & Xw)
{
	int n = Xw.size();
	point_t temp;
	points_t f;

	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2 * n, 3);
	Eigen::MatrixXd B1 = Eigen::MatrixXd::Zero(2 * n, 10);
	for (int j = 0; j < n; j++)
	{
		double uj = xn[j][0];
		double vj = xn[j][1];
		point_t pj = Xw[j];

		Eigen::MatrixXd Aj(2, 3);
		Aj << 1, 0, -uj, 0, 1, -vj;
		Eigen::MatrixXd B1j = MatrixCs(pj, uj, vj);

		A.row(2 * j) = Aj.row(0);
		A.row(2 * j + 1) = Aj.row(1);

		B1.row(2 * j) = B1j.row(0);
		B1.row(2 * j + 1) = B1j.row(1);
	}
	Eigen::MatrixXd A_temp = (A.transpose()*A).inverse()*A.transpose();
	Eigen::MatrixXd C1 = A_temp*B1;
	Eigen::MatrixXd E = A*C1 - B1;
	Eigen::MatrixXd G1 = E.transpose()*E;


	//calculate the normalized direction vector of image points;
	for (int i = 0; i < n;i++)
	{
		temp[0] = xn[i][0];
		temp[1] = xn[i][1];
		temp[2] = 1;
		temp = temp / temp.norm();
		f.push_back(temp);	
	}

	//select two points W1 and W2 whose projection length is the longest in image plane;
	int i1 = 0;
	int i2 = 1;
	double limn = f[i1][0] * f[i2][0] + f[i1][1] * f[i2][1] + f[i1][2] * f[i2][2];
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n;j++)
		{
			double L = f[i][0] * f[j][0] + f[i][1] * f[j][1] + f[i][2] * f[j][2];
			if (L<limn)
			{
				i1 = i;
				i2 = j;
				limn = L;
			}
		}
	}

	//create an intermediate frame;
	point_t P1 = Xw[i1];
	point_t P2 = Xw[i2];
	point_t P0 = (P1 + P2) / 2;
	point_t x = P2 - P0;
	x = x / x.norm();

	point_t temp1(0, 1, 0);
	point_t temp2(0, 0, 1);
	point_t y;
	point_t z;
	if (abs(temp1.transpose()*x)<abs(temp2.transpose()*x))
	{
		z = x.cross(temp1);
		z = z / z.norm();
		y = z.cross(x);
		y = y / y.norm();
	}
	else
	{
		y = x.cross(temp2);
		y = y / y.norm();
		z = x.cross(y);
		z = z / z.norm();
	}	
	rotation_t Ro;
	Ro.col(0) = x;
	Ro.col(1) = y;
	Ro.col(2) = z;

	// transform the target points from world frame into the intermediate frame by using Ro;
	points_t XX;
	for (int i = 0; i < n;i++)
	{
		point_t temp3=Ro.transpose()*(Xw[i] - P0);
		XX.push_back(temp3);
	}

	//Divide the n-point set into (n-2) 3-point subsets and setting up the P3P equations;

	point_t v1 = f[i1];
	point_t v2 = f[i2];
	double cg1 = v1.transpose()*v2;
	double sg1 = sqrt(1 - cg1*cg1);
	double D1 = (XX[i1] - XX[i2]).norm();

	Eigen::MatrixXd D4 = Eigen::MatrixXd::Zero(n-2,5);
	int j = 0;
	for (int i = 0; i < n;i++)
	{
		if (i==i1||i==i2)
		{
			continue;
		}
		
		point_t vi = f[i];
		double cg2 = v1.transpose()*vi;
		double cg3 = v2.transpose()*vi;
		double sg2 = sqrt(1 - cg2*cg2);
		double D2 = (XX[i1] - XX[i]).norm();
		double D3 = (XX[i] - XX[i2]).norm();

		//get the coefficient of the P3P equation from each subset;
		D4.row(j) = getp3p(cg1, cg2, cg3, sg1, sg2, D1, D2, D3);

		j = j + 1;
	}

	//get the 7th order polynomial, the deviation of the cos function;
	Eigen::MatrixXd D7 = Eigen::MatrixXd::Zero(1, 8);
	for (int i = 0; i < n - 2;i++)
	{
		D7 = D7 + getpoly7(D4.row(i));
	}

	std::vector<double> realRoots = gv::math::on_rootss(D7);

	//calculate the camera pose from each local minimal;
	// get complete pose;
	transformations_t pose;

	std::vector<double> error;
	for (unsigned i = 0; i < realRoots.size();i++)
	{
		double t2 = realRoots[i];
		double d2 = cg1 + t2;
		x = v2*d2 - v1;
		x = x / x.norm();

		if (abs(temp1.transpose()*x)<abs(temp2.transpose()*x))
		{
			z = x.cross(temp1);
			z = z / z.norm();
			y = z.cross(x);
			y = y / y.norm();
		}
		else
		{
			y = temp2.cross(x);
			y = y / y.norm();
			z = x.cross(y);
			z = z / z.norm();
		}

		rotation_t Rx;
		Rx.col(0) = x;
		Rx.col(1) = y;
		Rx.col(2) = z;

		Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2*n, 3);
		for (int j = 0; j < n;j++)
		{
			double uj = xn[j][0];
			double vj = xn[j][1];
			point_t Pj = XX[j];

			Eigen::MatrixXd Bj = fucntionB(Rx,Pj,uj,vj);

			B.row(2 * j) = Bj.row(0);
			B.row(2 * j + 1) = Bj.row(1);
		}
		Eigen::MatrixXd C = A_temp*B;
		Eigen::MatrixXd E1 = A*C - B;
		Eigen::MatrixXd G = E1.transpose()*E1;

		double g11 = G(0, 0); double g12 = G(0, 1); double g13 = G(0, 2);
		double g22 = G(1, 1); double g23 = G(1, 2);

		Eigen::Matrix<double, 5, 1> factors;
		factors(0, 0) = 4 * pow(g12, 2) + pow(g22, 2) + pow(g11, 2) - 2 * g11*g22;
		factors(1, 0) = 4 * g12*g23 + 2 * g11*g13 - 2 * g13*g22;
		factors(2, 0) = pow(g23, 2) + 2 * g11*g22 + pow(g13, 2) - 4 * pow(g12, 2) - pow(g11, 2) - pow(g22, 2);
		factors(3, 0) = 2 * g13*g22 - 2 * g11*g13 - 2 * g12*g23;
		factors(4, 0) = g12*g12 - g13*g13;

		std::vector<double> realRoots1 = gv::math::o4_roots(factors);

		for (unsigned j = 0; j < realRoots1.size(); j++)
		{
			double cj = realRoots1[j];
			double sj = (2 * g12*pow(cj, 2) + g23*cj - g12) / ((g11 - g22)*cj + g13);
			point_t ss(cj,sj,1);
			point_t ts = C*ss;

			double r1 = Rx(0, 0); double r2 = Rx(0, 1); double r3 = Rx(0, 2);
			double r4 = Rx(1, 0); double r5 = Rx(1, 1); double r6 = Rx(1, 2);
			double r7 = Rx(2, 0); double r8 = Rx(2, 1); double r9 = Rx(2, 2);
			
			//calculate the camera pose by 3d alignment;
			points_t XXc;
			for (int k = 0; k < n;k++)
			{
				double xi = XX[k][0];
				double yi = XX[k][1];
				double zi = XX[k][2];
				double x_temp = r1*xi + (r2*cj + r3*sj)*yi + (-r2*sj + r3*cj)*zi + ts[0];
				double y_temp = r4*xi + (r5*cj + r6*sj)*yi + (-r5*sj + r6*cj)*zi + ts[1];
				double z_temp = r7*xi + (r8*cj + r9*sj)*yi + (-r8*sj + r9*cj)*zi + ts[2];
				point_t XXcs_temp(x_temp,y_temp,z_temp);
				XXc.push_back(f[k] * XXcs_temp.norm());
			}
			transformation_t solution = math::arun_complete(Xw, XXc);

			rotation_t R_estimation = solution.block<3, 3>(0, 0).transpose();
			translation_t t_estimation = -solution.block<3, 3>(0, 0).transpose()*solution.col(3);
			double er = 0;
			for (int k = 0; k < n;k++)
			{
				point_t Xc = R_estimation*Xw[k] + t_estimation;
				Image_point xc(Xc[0] / Xc[2], Xc[1] / Xc[2]);
				er = er+(xc - xn[k]).norm();
			}
			error.push_back(er / n);
			transformation_t solution1;
			solution1.block<3, 3>(0, 0) = R_estimation;
			solution1.col(3) = t_estimation;
			pose.push_back(solution1);

		}
	}//solution end;

	//determine the camera pose with the smallest projection error;
	double minr = error[0];
	int index = 0;
	for (unsigned j = 0; j < error.size();j++)
	{
		if (error[j]<minr)
		{
			minr = error[j];
			index = j;
		}
	}

	transformation_t pose1;
	pose1 = pose[index];
	
	//To further improve the accuracy by using a single Gauss_Newton method;
	//convert rotation matrix to quaternion;
	//Eigen::Vector4d sol2 = matrix2quaternion(pose1.block<3, 3>(0, 0));
	//Eigen::Vector3d sol3(sol2[1] / sol2[0], sol2[2] / sol2[0], sol2[3] / sol2[0]);
	
	Eigen::Vector3d sol2 = Cayley(pose1.block<3, 3>(0, 0));

	Eigen::Vector3d solution=RefineGaussNewton(sol2,E,G1);
	double s1 = solution[0];
	double s2 = solution[1];
	double s3 = solution[2];
	
	rotation_t Rw;
	Rw << 1 + s1*s1 - s2 *s2 - s3 *s3, 2 * s1*s2 - 2 * s3, 2 * s2 + 2 * s1*s3,
		  2 * s3 + 2 * s1*s2, 1 - s1 *s1 + s2 *s2 - s3 *s3, 2 * s2*s3 - 2 * s1,
		  2 * s1*s3 - 2 * s2, 2 * s1 + 2 * s2*s3, 1 - s1 *s1 - s2 *s2 + s3 *s3;
	
	Eigen::VectorXd w(10);
	w << 1, s1, s2, s3, s1*s1, s1*s2, s1*s3, s2*s2, s2*s3, s3 *s3;
	double factor = 1 / (1 + s1*s1 + s2*s2 + s3*s3);

	transformation_t pose3;
	pose3.block<3, 3>(0, 0) = Rw*factor;
	pose3.col(3) = C1*w*factor;

	return pose3;

}//SRPnP end;

Eigen::MatrixXd
gv::SRPnP::getp3p(double l1, double l2, double A5, double C1, double C2, double D1, double D2, double D3)
{
	double A1 = pow(D2 / D1,2);
	double A2 = A1*C1*C1 - C2*C2;
	double A3 = l2*A5 - l1;
	double A4 = l1*A5 - l2;
	double A6 = (D3*D3 - D1*D1 - D2*D2) / (2 * D1*D1);
	double A7 = 1 - l1*l1 - l2*l2 + l1*l2*A5 + A6*C1*C1;

	Eigen::MatrixXd B(1, 5);
	B << A6*A6 - A1*A5*A5, 2 * (A3*A6 - A1*A4*A5), A3*A3 + 2 * A6*A7 - A1*A4*A4 - A2*A5*A5,
		2 * (A3*A7 - A2*A4*A5), A7*A7 - A2*A4*A4;

	return B;
}

Eigen::MatrixXd
gv::SRPnP::getpoly7(Eigen::MatrixXd F)
{
	Eigen::MatrixXd F7(1, 8);
	F7 << 4 * F(0)*F(0), 7 * F(1)*F(0), 6 * F(2)*F(0) + 3 * F(1)*F(1), 5 * F(3)*F(0) + 5 * F(2)*F(1),
		  4 * F(4)*F(0) + 4 * F(3)*F(1) + 2 * F(2)*F(2), 3 * F(4)*F(1) + 3 * F(3)*F(2), 2 * F(4)*F(2) + F(3)*F(3), F(4)*F(3);

	return F7;
}

Eigen::MatrixXd
gv::SRPnP::fucntionB(rotation_t R, point_t P, double ui, double vi)
{
	double Xi = P[0];	double Yi = P[1];	double Zi=P[2];
	double r1 = R(0, 0); double r2 = R(0, 1); double r3 = R(0, 2);
	double r4 = R(1, 0); double r5 = R(1, 1); double r6 = R(1, 2);
	double r7 = R(2, 0); double r8 = R(2, 1); double r9 = R(2, 2);

	Eigen::MatrixXd B(2, 3);
	B << Yi*r8*ui + Zi*r9*ui - Yi*r2 - Zi*r3, Yi*r9*ui - Zi*r8*ui - Yi*r3 + Zi*r2, Xi*r7*ui - Xi*r1,
		 Yi*r8*vi + Zi*r9*vi - Yi*r5 - Zi*r6, Yi*r9*vi - Zi*r8*vi - Yi*r6 + Zi*r5, Xi*r7*vi - Xi*r4;

	return B;
}

Eigen::MatrixXd
gv::SRPnP::MatrixCs(point_t P, double ui, double vi)
{
	double X = P[0];	double Y = P[1];	double Z = P[2];

	Eigen::MatrixXd C(2, 10);
	C << -X + Z*ui, 2 * Y*ui, -2 * Z - 2 * X*ui, 2 * Y, -X - Z*ui, -2 * Y, -2 * Z + 2 * X*ui, X - Z*ui, 2 * Y*ui, Z*ui + X,
		 -Y + Z*vi, 2 * Y*vi + 2 * Z, -2 * X*vi, -2 * X, Y - Z*vi, -2 * X, 2 * X*vi, -Y - Z*vi, -2 * Z + 2 * Y*vi, Z*vi + Y;

	return C;
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

//Convert rotation matrix to quaternion
Eigen::Vector4d
gv::SRPnP::matrix2quaternion(rotation_t R)
{
	Eigen::MatrixXd I;
	I.setIdentity(3,3);

	Eigen::EigenSolver<Eigen::Matrix3d> es(R - I);
	//Eigen::Matrix3d D = es.pseudoEigenvalueMatrix();
	Eigen::Matrix3d V = es.pseudoEigenvectors();

	std::vector<double> d;
	for (int i = 0; i < es.eigenvalues().size();i++)
	{
		d.push_back(es.eigenvalues().row(i).norm());
	}
	double mind = d[0];
	int index1 = 0;
	for (unsigned i = 0; i < d.size();i++)
	{
		if (d[i]<mind)
		{
			mind = d[i];
			index1 = i;
		}
	}
	point_t axis = V.col(index1);

	double twocostheta = R.trace() - 1;
	point_t twosinthetav(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
	double twosintheta = axis.transpose()*twosinthetav;
	double theta = atan2(twosintheta, twocostheta);

	Eigen::Vector4d Q;
	Q(0) = cos(theta / 2);
	Q.segment(1, 3) = axis*sin(theta/2);


	return Q;
}

//Convert rotation matrix to quaternion
Eigen::Vector3d
gv::SRPnP::Cayley(rotation_t R)
{
	double q1, q2, q3, q4;
	rotation_t A = R.transpose();
	q4 = sqrt(1 + A(0, 0) + A(1, 1) + A(2, 2)) / 2;

	if (q4 > 0.01)
	{
		q1 = (A(2, 1) - A(1, 2)) *0.25 / q4;
		q2 = (A(0, 2) - A(2, 0)) *0.25 / q4;
		q3 = (A(1, 0) - A(0, 1)) *0.25 / q4;
	}
	else
	{
		q1 = sqrt(1 + A(0, 0) - A(1, 1) - A(2, 2)) / 2;
		q2 = (A(0, 1) + A(1, 0))*0.25 / q1;
		q3 = (A(0, 2) + A(2, 0))*0.25 / q1;
		q4 = (A(2, 1) - A(1, 2))*0.25 / q1;
	}

	Eigen::Vector3d Q;
	Q(0) = -q1 / q4;
	Q(1) = -q2 / q4;
	Q(2) = -q3 / q4;

	return Q;
}

//Refine the solution by using one step Gauss Newton method;
//By Ping Wang
Eigen::Vector3d
gv::SRPnP::RefineGaussNewton(Eigen::Vector3d solution, Eigen::MatrixXd E, Eigen::MatrixXd G1)
{
	double s1 = solution[0];
	double s2 = solution[1];
	double s3 = solution[2];

	Eigen::VectorXd w(10);
	w << 1, s1, s2, s3, s1*s1, s1*s2, s1*s3, s2*s2, s2*s3, s3 *s3;
	double obj_pre = w.transpose()*G1*w;
	//Jacobian matrix;
	Eigen::MatrixXd Jac(3, 10);
	Jac << 0, 1, 0, 0, 2 * s1, s2, s3, 0, 0, 0,
		   0, 0, 1, 0, 0, s1, 0, 2 * s2, s3, 0,
		   0, 0, 0, 1, 0, 0, s1, 0, s2, 2 * s3;
	Eigen::MatrixXd Fk = E*w;
	Eigen::MatrixXd Jk = E*Jac.transpose();
	Eigen::Vector3d solution_temp = solution;
	//increment
	Eigen::Vector3d dk = ((Jk.transpose()*Jk).inverse()*Jk.transpose()*Fk);
	//update parameter
	solution_temp = solution_temp - dk;
	s1 = solution_temp[0];
	s2 = solution_temp[1];
	s3 = solution_temp[2];
	w << 1, s1, s2, s3, s1*s1, s1*s2, s1*s3, s2*s2, s2*s3, s3 *s3;
	//evaluate the error of objection
	double obj_cur = w.transpose()*G1*w;
	if (obj_cur<obj_pre)
	{
		solution = solution_temp;
	}

	return solution;
}


