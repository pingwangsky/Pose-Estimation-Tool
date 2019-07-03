/******************************************************************************
* This programe is implemented in Visual Studio 2013.                        *
* Author :  Ping Wang                                                        *
* Contact:  pingwangsky@gmail.com                                            *
* License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
* Address:  Lanzhou University of Technology                                 *
* My site:  https://sites.google.com/view/ping-wang-homepage                 *
******************************************************************************/

/**
 * \file types.hpp
 * \brief A collection of variables used in geometric vision for the
 *        computation of calibrated absolute and relative pose.
 */

#ifndef GV_TYPES_HPP_
#define GV_TYPES_HPP_

#include <stdlib.h>
#include <vector>
#include <Eigen3\Eigen\Eigen>
#include <Eigen3\Eigen\Dense>
#include <Eigen3\Eigen\src\Core\util\DisableStupidWarnings.h>


/**
 * \brief The namespace of this library.
 */
namespace gv
{

/** A 3-vector of unit length used to describe landmark observations/bearings
 *  in camera frames (always expressed in camera frames)
 */
typedef Eigen::Vector3d
    bearingVector_t;

/** An array of bearing-vectors */
typedef std::vector<bearingVector_t, Eigen::aligned_allocator<bearingVector_t> >
    bearingVectors_t;

/** A 3-vector describing a translation/camera position */
typedef Eigen::Vector3d
    translation_t;

/** An array of translations */
typedef std::vector<translation_t, Eigen::aligned_allocator<translation_t> >
    translations_t;

/** A rotation matrix */
typedef Eigen::Matrix3d
    rotation_t;

/** An array of rotation matrices as returned by fivept_kneip [7] */
//需要16字节对齐; 
typedef std::vector<rotation_t, Eigen::aligned_allocator<rotation_t> >
    rotations_t;

/** A 3x4 transformation matrix containing rotation \f$ \mathbf{R} \f$ and
 *  translation \f$ \mathbf{t} \f$ as follows:
 *  \f$ \left( \begin{array}{cc} \mathbf{R} & \mathbf{t} \end{array} \right) \f$
 */
typedef Eigen::Matrix<double,3,4>
    transformation_t;

/** An array of transformations */
typedef std::vector<transformation_t, Eigen::aligned_allocator<transformation_t> >
    transformations_t;

/** A 3-vector containing the cayley parameters of a rotation matrix */
typedef Eigen::Vector3d
    cayley_t;

/** A 4-vector containing the quaternion parameters of rotation matrix */
typedef Eigen::Vector4d
    quaternion_t;

/** Essential matrix \f$ \mathbf{E} \f$ between two viewpoints:
 *
 *  \f$ \mathbf{E} = \f$ skew(\f$\mathbf{t}\f$) \f$ \mathbf{R} \f$,
 *
 *  where \f$ \mathbf{t} \f$ describes the position of viewpoint 2 seen from
 *  viewpoint 1, and \f$\mathbf{R}\f$ describes the rotation from viewpoint 2
 *  to viewpoint 1.
 */
typedef Eigen::Matrix3d
    essential_t;

/** An array of essential matrices */
typedef std::vector<essential_t, Eigen::aligned_allocator<essential_t> >
    essentials_t;

/** An essential matrix with complex entires (as returned from
 *  fivept_stewenius [5])
 */
typedef Eigen::Matrix3cd
    complexEssential_t;

/** An array of complex-type essential matrices */
typedef std::vector< complexEssential_t, Eigen::aligned_allocator< complexEssential_t> >
    complexEssentials_t;

/** A 3-vector describing a point in 3D-space */
typedef Eigen::Vector2d
    Image_point;

/** An array of 3D-points */
typedef std::vector<Image_point, Eigen::aligned_allocator<Image_point> >
    Image_points;

/** A 3-vector describing a point in 3D-space */
typedef Eigen::Vector3d
    point_t;

/** An array of 3D-points */
typedef std::vector<point_t, Eigen::aligned_allocator<point_t> >
    points_t;

/** A 3-vector containing the Eigenvalues of matrix \f$ \mathbf{M} \f$ in the
 *  eigensolver-algorithm (described in [11])
 */
typedef Eigen::Vector3d
    eigenvalues_t;

/** A 3x3 matrix containing the eigenvectors of matrix \f$ \mathbf{M} \f$ in the
 *  eigensolver-algorithm (described in [11])
 */
typedef Eigen::Matrix3d
    eigenvectors_t;

/*定义实数*/
typedef std::vector<double, Eigen::aligned_allocator<double> > Coeff;

/** 定义复数
*/
typedef std::complex<double> complex_c;

/** An array of transformations */
typedef std::vector<complex_c, Eigen::aligned_allocator<complex_c> > complex_C;

/** EigensolverOutput holds the output-parameters of the eigensolver-algorithm
 *  (described in [11])
 */
typedef struct EigensolverOutput
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** Position of viewpoint 2 seen from viewpoint 1 (unscaled) */
  translation_t   translation;
  /** Rotation from viewpoint 2 back to viewpoint 1 */
  rotation_t      rotation;
  /** The eigenvalues of matrix \f$ \mathbf{M} \f$ */
  eigenvalues_t   eigenvalues;
  /** The eigenvectors of matrix matrix \f$ \mathbf{M} \f$ */
  eigenvectors_t  eigenvectors;
} eigensolverOutput_t;

/** GeOutput holds the output-parameters of ge
 */
typedef struct GeOutput
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  /** Homogeneous position of viewpoint 2 seen from viewpoint 1 */
  Eigen::Vector4d   translation;
  /** Rotation from viewpoint 2 back to viewpoint 1 */
  rotation_t        rotation;
  /** The eigenvalues of matrix \f$ \mathbf{G} \f$ */
  Eigen::Vector4d   eigenvalues;
  /** The eigenvectors of matrix matrix \f$ \mathbf{G} \f$ */
  Eigen::Matrix4d   eigenvectors;
} geOutput_t;

}

#endif /* GV_TYPES_HPP_ */
