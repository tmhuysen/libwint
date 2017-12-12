#ifndef LIBWRP_TRANSFORMATIONS_HPP
#define LIBWRP_TRANSFORMATIONS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


namespace libwrp {


/** Given:
 *      - a matrix h, which contains one-electron integrals in some basis B
 *      - a transformation matrix which represents the transformation of B into B'
 *
 * return the matrix of one-electron integrals h' in the new basis B'.
 *
 *
 *
 *  Note that the transformation matrix T is used as
 *
 *      B' = B T ,
 *
 *  where the basis vectors are collected as elements of a row vector.
 */
Eigen::MatrixXd transform_one_electron_integrals(Eigen::MatrixXd& h, Eigen::MatrixXd& T);


/** Given:
 *      - a rank-four tensor g, which contains the two-electron integrals in some basis B
 *      - a transformation matrix which represents the transformation of B into B'
 *
 * return the tensor g' in the new basis B'.
 *
 *
 *
 *  Note that the transformation matrix T is used as
 *
 *      B' = B T ,
 *
 *  where the basis vectors are collected as elements of a row vector.
 */
Eigen::Tensor<double, 4> transform_two_electron_integrals(Eigen::Tensor<double, 4>& g, Eigen::MatrixXd& T);


/** Given:
 *      - a matrix representation in an AO basis (f_AO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the SO basis
 */
Eigen::MatrixXd transform_AO_to_SO(Eigen::MatrixXd& f_AO, Eigen::MatrixXd& C);


/** Given:
 *      - a matrix representation in an SO basis (f_SO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the AO basis
 */
Eigen::MatrixXd transform_SO_to_AO(Eigen::MatrixXd& f_SO, Eigen::MatrixXd& C);


/** Given:
 *      - a rank-four tensor of two-electron integrals in an AO basis
 *      - an SO coefficient matrix (every column represents a spatial orbital)
 *
 *  transform and return the two-electron integrals in the SO basis
 */
Eigen::Tensor<double, 4> transform_AO_to_SO(Eigen::Tensor<double, 4>& g_AO, Eigen::MatrixXd& C);


/** Given a unitary matrix U that transforms a basis B into B', return the one-electron integrals in the rotated basis
 *
 * Note that the basis transformation is explicitly written as (B' = B U)
 */
Eigen::MatrixXd rotate_integrals(Eigen::MatrixXd& h, Eigen::MatrixXd& U);


/** Given a unitary matrix U that transforms a basis B into B', return the two-electron integrals in the rotated basis
 *
 * Note that the basis transformation is explicitly written as (B' = B U)
 */
Eigen::Tensor<double, 4> rotate_integrals(Eigen::Tensor<double, 4>& g, Eigen::MatrixXd& U);


/** Give the M-dimensional Jacobi rotation matrix (with an angle theta) for the orbitals P and Q
 *
 * Note that we work with the (cos, sin, -sin, cos) definition
 */
Eigen::MatrixXd jacobi_rotation_matrix(size_t P, size_t Q, double theta, size_t M);


}  // namespace libwrp

#endif // LIBWRP_TRANSFORMATIONS_HPP
