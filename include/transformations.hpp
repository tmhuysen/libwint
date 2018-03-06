#ifndef LIBWINT_TRANSFORMATIONS_HPP
#define LIBWINT_TRANSFORMATIONS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>



namespace libwint {
namespace transformations {


/*
 *  GENERAL TRANSFORMATIONS
 */
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
Eigen::MatrixXd transformOneElectronIntegrals(const Eigen::MatrixXd& h, const Eigen::MatrixXd& T);

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
Eigen::Tensor<double, 4> transformTwoElectronIntegrals(const Eigen::Tensor<double, 4>& g, const Eigen::MatrixXd& T);



/*
 *  AO AND SO CONVERSION WRAPPERS
 */

/** Given:
 *      - a matrix representation in an AO basis (f_AO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the SO basis
 */
Eigen::MatrixXd transform_AO_to_SO(const Eigen::MatrixXd& f_AO, const Eigen::MatrixXd& C);

/** Given:
 *      - a matrix representation in an SO basis (f_SO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the AO basis
 */
Eigen::MatrixXd transform_SO_to_AO(const Eigen::MatrixXd& f_SO, const Eigen::MatrixXd& C);

/** Given:
 *      - a rank-four tensor of two-electron integrals in an AO basis
 *      - an SO coefficient matrix (every column represents a spatial orbital)
 *
 *  transform and return the two-electron integrals in the SO basis
 */
Eigen::Tensor<double, 4> transform_AO_to_SO(const Eigen::Tensor<double, 4>& g_AO, const Eigen::MatrixXd& C);


/*
 *  JACOBI ROTATIONS AND WRAPPERS
 */

/** Give the M-dimensional Jacobi rotation matrix for the orbitals p and q (p < q) and a given @param theta.
 *
 * M is the actual dimension of the matrix that is returned
 * @param p and @param q represent the rows and columns, i.e. they start at 0
 *
 * Note that we work with the (cos, sin, -sin, cos) definition
 */
Eigen::MatrixXd jacobiRotationMatrix(size_t p, size_t q, double theta, size_t M);

/**
 *  Using a Jacobi rotation with angle @param: theta (in radians) of the orbitals p and q, return the transformed one-electron integrals.
 *  This function is implemented using Eigen's Jacobi module.
 */
Eigen::MatrixXd rotateOneElectronIntegralsJacobi(const Eigen::MatrixXd& h, size_t p, size_t q, double theta);

/**
 *  Using a Jacobi rotation with angle @param: theta (in radians) of the orbitals p and q, return the transformed two-electron integrals.
 *  While waiting for an analogous Eigen::Tensor Jacobi module, this function is just a wrapper around transformTwoElectronIntegrals using a Jacobi rotation matrix.
 */
Eigen::Tensor<double, 4> rotateTwoElectronIntegralsJacobi(const Eigen::Tensor<double, 4>& g, size_t p, size_t q, double theta);



}  // namespace transformations
}  // namespace libwint

#endif // LIBWINT_TRANSFORMATIONS_HPP
