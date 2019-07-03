/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/

/**
 * \file roots.hpp
 * \brief Closed-form solutions for computing the roots of a polynomial.
 */

#ifndef GV_ROOTS_HPP_
#define GV_ROOTS_HPP_

#include <stdlib.h>
#include <vector>
#include <Eigen3/Eigen/Eigen>
#include <Eigen3/Eigen/src/Core/util/DisableStupidWarnings.h>
#include "gv\types.hpp"
/**
 * \brief The namespace of this library.
 */
namespace gv
{
/**
 * \brief The namespace of the math tools.
 */
namespace math
{
/**
 * \brief The real roots of a second-order polynomial.
 * 
 * \param[in] p The polynomial coefficients (poly = p[0,0]*x^2 + p[1,0]*x^2+p[2,0] ...).
 * \return The roots of the polynomial (only real ones).
 */
std::vector<double> o2_roots( const Eigen::MatrixXd & p );

/**
 * \brief The roots of a third-order polynomial.
 * 
 * \param[in] p The polynomial coefficients (poly = p[0]*x^3 + p[1]*x^2 ...).
 * \return The roots of the polynomial (only real ones).
 */
std::vector<double> o3_roots( const std::vector<double> & p );

/**
 * \brief Ferrari's method for computing the roots of a fourth order polynomial.
 *
 * \param[in] p The polynomial coefficients (poly = p(0,0)*x^4 + p(1,0)*x^3 ...).
 * \return The roots of the polynomial (only real ones).
 */
std::vector<double> o4_roots( const Eigen::MatrixXd & p );

/**
 * \brief Ferrari's method for computing the roots of a fourth order polynomial.
 *        With a different interface.
 *
 * \param[in] p The polynomial coefficients (poly = p[0]*x^4 + p[1]*x^3 ...).
 * \return The roots of the polynomial (only real ones).
 */
std::vector<double> o4_roots( const std::vector<double> & p );

//输入是动态向量形式
gv::complex_C const on_roots(const gv::Coeff sol);


//
std::vector<double> on_rootss(const Eigen::MatrixXd & sol);

}
}

#endif /* GV_ROOTS_HPP_ */
