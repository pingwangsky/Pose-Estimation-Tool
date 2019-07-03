
/**
 * \file arun.hpp
 * \brief Arun's method for computing the rotation between two point sets.
 */

#ifndef GV_ARUN_HPP_
#define GV_ARUN_HPP_

#include <stdlib.h>
#include <Eigen3/Eigen/Eigen>
#include <Eigen3/Eigen/src/Core/util/DisableStupidWarnings.h>
#include <gv/types.hpp>

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
 * \brief Arun's method for computing the rotation between two point sets.
 *        Core function [13].
 *
 * \param[in] Hcross The summation over the exterior products between the
 *            normalized points.
 * \return The rotation matrix that aligns the points.
 */
rotation_t arun( const Eigen::MatrixXd & Hcross );

/**
 * \brief Arun's method for complete point cloud alignment [13]. The method
 *        actually does the same than threept_arun, but has a different
 *        interface.
 *
 * \param[in] p1 The points expressed in the first frame.
 * \param[in] p2 The points expressed in the second frame.
 * \return The Transformation from frame 2 to frame 1 (
 *         \f$ \mathbf{T} = \left(\begin{array}{cc} \mathbf{R} & \mathbf{t} \end{array}\right) \f$,
 *         with \f$ \mathbf{t} \f$ being the position of frame 2 seen from
 *         frame 1, and \f$ \mathbf{R} \f$ being the rotation from
 *         frame 2 to frame 1).
 */
transformation_t arun_complete( const points_t & p1, const points_t & p2 );

}
}

#endif /* OPENGV_ARUN_HPP_ */
