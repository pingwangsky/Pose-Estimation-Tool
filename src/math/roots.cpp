/******************************************************************************
* This programe is implemented in Visual Studio 2015.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/

#include <gv/math/roots.hpp>
#include <complex>

//calculate the real roots of second-order ploynomial;
std::vector<double>
gv::math::o2_roots( const Eigen::MatrixXd & p )
{
  double A = p(0,0);
  double B = p(1,0);
  double C = p(2,0);
  double r1;
  double r2;

  double delt=pow(B,2)-4*A*C;
  if(delt<=0)
  {
	  r1=-B/(2*A);
	  r2=-B/(2*A);
  }
  else if(delt>0)
  { 
	  r1=(-B+sqrt(delt))/(2*A);
	  r2=(-B-sqrt(delt))/(2*A);
  }

  std::vector<double> realRoots;
  realRoots.push_back(r1);
  realRoots.push_back(r2);

  return realRoots;
}

std::vector<double>
gv::math::o3_roots( const std::vector<double> & p )
{
  const double & a = p[0];
  const double & b = p[1];
  const double & c = p[2];
  const double & d = p[3];
  
  double theta_0 = b*b - 3.0*a*c;
  double theta_1 = 2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d;
  double term = theta_1 * theta_1 - 4.0 * theta_0 * theta_0 * theta_0;
  
  std::complex<double> u1( 1.0, 0.0);
  std::complex<double> u2(-0.5, 0.5*sqrt(3.0));
  std::complex<double> u3(-0.5,-0.5*sqrt(3.0));
  std::complex<double> C;
  
  if( term >= 0.0 )
  {
    double C3 = 0.5 * (theta_1 + sqrt(term));
    
    if( C3 < 0.0 )
      C = std::complex<double>(-pow(-C3,(1.0/3.0)),0.0);
    else
      C = std::complex<double>( pow( C3,(1.0/3.0)),0.0);
  }
  else
  {
    std::complex<double> C3( 0.5*theta_1, 0.5*sqrt(-term) );
    
    //take the third root of this complex number
    double r3 = sqrt(pow(C3.real(),2.0)+pow(C3.imag(),2.0));
    double a3 = atan(C3.imag() / C3.real());
    if( C3.real() < 0 )
      a3 += M_PI;
    
    double r = pow(r3,(1.0/3.0));
    C = std::complex<double>(r*cos(a3/3.0),r*sin(a3/3.0));
  }
  
  std::complex<double> r1 = -(1.0/(3.0*a)) * (b + u1*C + theta_0 / (u1*C) );
  std::complex<double> r2 = -(1.0/(3.0*a)) * (b + u2*C + theta_0 / (u2*C) );
  std::complex<double> r3 = -(1.0/(3.0*a)) * (b + u3*C + theta_0 / (u3*C) );
  
  std::vector<double> roots;
  roots.push_back(r1.real());
  roots.push_back(r2.real());
  roots.push_back(r3.real());
  return roots;
}

std::vector<double>
gv::math::o4_roots( const Eigen::MatrixXd & p )
{
  double A = p(0,0);
  double B = p(1,0);
  double C = p(2,0);
  double D = p(3,0);
  double E = p(4,0);

  double A_pw2 = A*A;
  double B_pw2 = B*B;
  double A_pw3 = A_pw2*A;
  double B_pw3 = B_pw2*B;
  double A_pw4 = A_pw3*A;
  double B_pw4 = B_pw3*B;

  double alpha = -3*B_pw2/(8*A_pw2)+C/A;
  double beta = B_pw3/(8*A_pw3)-B*C/(2*A_pw2)+D/A;
  double gamma = -3*B_pw4/(256*A_pw4)+B_pw2*C/(16*A_pw3)-B*D/(4*A_pw2)+E/A;

  double alpha_pw2 = alpha*alpha;
  double alpha_pw3 = alpha_pw2*alpha;

  std::complex<double> P (-alpha_pw2/12-gamma,0);
  std::complex<double> Q (-alpha_pw3/108+alpha*gamma/3-pow(beta,2)/8,0);
  std::complex<double> R = -Q/2.0+sqrt(pow(Q,2.0)/4.0+pow(P,3.0)/27.0);

  std::complex<double> U = pow(R,(1.0/3.0));
  std::complex<double> y;

  if (U.real() == 0)
    y = -5.0*alpha/6.0-pow(Q,(1.0/3.0));
  else
    y = -5.0*alpha/6.0-P/(3.0*U)+U;

  std::complex<double> w = sqrt(alpha+2.0*y);

  std::vector<double> realRoots;
  std::complex<double> temp;
  temp = -B/(4.0*A) + 0.5*(w+sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
  if (abs(temp.imag()) <0.001)
  {
	  realRoots.push_back(temp.real());
  }
  //realRoots.push_back(temp.real());
  temp = -B/(4.0*A) + 0.5*(w-sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
  if (abs(temp.imag()) <0.001)
  {
	  realRoots.push_back(temp.real());
  }
  //realRoots.push_back(temp.real());
  temp = -B/(4.0*A) + 0.5*(-w+sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
  if (abs(temp.imag()) <0.001)
  {
	  realRoots.push_back(temp.real());
  }
  //realRoots.push_back(temp.real());
  temp = -B/(4.0*A) + 0.5*(-w-sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
  if (abs(temp.imag()) <0.001)
  {
	  realRoots.push_back(temp.real());
  }
  //realRoots.push_back(temp.real());

  return realRoots;
}


//calculate the real roots of n-order ploynomial;
gv::complex_C const
gv::math::on_roots(const gv::Coeff sol)
{
	int n = sol.size();
	Eigen::VectorXd vec(n - 1);
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n-1, n-1);

	for (int i = 0; i <n-1;i++)
	{
		vec[i] = sol[i + 1] / sol[0];
	}
	M.row(0) = -vec;
	for (unsigned int i = 0; i < n - 2;i++)
	{
		M(i + 1, i) = 1;
	}

	gv::complex_C solution;
	Eigen::MatrixXcd E = M.eigenvalues();
	
	for (int i = 0; i < E.size(); i++)
	{
		solution.push_back(E.row(i)[0]);
	}

	return solution;
}


std::vector<double>
gv::math::on_rootss(const Eigen::MatrixXd & sol)
{
	int n = sol.cols();
	//Eigen::VectorXd vec = sol.row(0).rightCols(n - 1) / sol(0, 0);
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n - 1, n - 1);

	Eigen::VectorXd vec(n - 1);
	for (int i = 0; i <n - 1; i++)
	{
		vec[i] = sol(0,i+1)/ sol(0,0);
	}

	M.row(0) = -vec;
	for (unsigned int i = 0; i < n - 2; i++)
	{
		M(i + 1, i) = 1;
	}
	Eigen::MatrixXcd E = M.eigenvalues();
	//double maxreal = E.real().cwiseAbs().maxCoeff();
	//select minima;
	double maxreal = abs(E.row(0)[0].real());
	for (unsigned int i = 0; i < E.size(); i++)
	{
		if (E.row(i)[0].real()>maxreal)
		{
			maxreal = abs(E.row(i)[0].real());
		}
	}
	std::vector<double> realRoots;
	realRoots.push_back(0);
	for (int i = 0; i < E.size(); i++)
	{
		if (abs(E.row(i)[0].imag()) / maxreal<0.001)
		{
			realRoots.push_back(E.row(i)[0].real());
		}
	}

	return realRoots;

}