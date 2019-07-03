
#include <gv/math/arun.hpp>

gv::rotation_t
gv::math::arun( const Eigen::MatrixXd & Hcross )
{
  //decompose matrix H to obtain rotation
  Eigen::JacobiSVD< Eigen::MatrixXd > SVDcross(
      Hcross,
      Eigen::ComputeFullU | Eigen::ComputeFullV );

  Eigen::Matrix3d V = SVDcross.matrixV();
  Eigen::Matrix3d U = SVDcross.matrixU();
  rotation_t R = V * U.transpose();

  //modify the result in case the rotation has determinant=-1
  if( R.determinant() < 0 )
  {
    Eigen::Matrix3d V_prime;
    V_prime.col(0) = V.col(0);
    V_prime.col(1) = V.col(1);
    V_prime.col(2) = -V.col(2);
    R = V_prime * U.transpose();
  }

  return R;
}

gv::transformation_t
gv::math::arun_complete(
    const points_t & p1,
    const points_t & p2 )
{
  assert(p1.size() == p2.size());

  //derive the centroid of the two point-clouds
  point_t pointsCenter1 = Eigen::Vector3d::Zero();
  point_t pointsCenter2 = Eigen::Vector3d::Zero();

  for( size_t i = 0; i < p1.size(); i++ )
  {
    pointsCenter1 += p1[i];
    pointsCenter2 += p2[i];
  }

  pointsCenter1 = pointsCenter1 / p1.size();
  pointsCenter2 = pointsCenter2 / p2.size();

  //compute the matrix H = sum(f'*f^{T})
  Eigen::MatrixXd Hcross(3,3);
  Hcross = Eigen::Matrix3d::Zero();

  for( size_t i = 0; i < p1.size(); i++ )
  {
    Eigen::Vector3d f = p1[i] - pointsCenter1;
    Eigen::Vector3d fprime = p2[i] - pointsCenter2;
    Hcross += fprime * f.transpose();
  }

  //decompose this matrix (SVD) to obtain rotation
  rotation_t rotation = arun(Hcross);
  translation_t translation = pointsCenter1 - rotation*pointsCenter2;
  transformation_t solution;
  solution.block<3,3>(0,0) = rotation;
  solution.col(3) = translation;

  return solution;
}
