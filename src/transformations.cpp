#include "transformations.hpp"

#include <iostream>

#include <Eigen/Jacobi>



namespace libwint {
namespace transformations {


/*
 *  GENERAL TRANSFORMATIONS
 */

/** Given:
 *      - a matrix h, which contains one-electron integrals in some basis B
 *      - a transformation matrix which represents the transformation of B into B'
 *
 * return the matrix of one-electron integrals h' in the new basis B'
 *
 *
 *
 *
 *  Note that the transformation matrix T is used as
 *
 *      B' = B T ,
 *
 *  where the basis vectors are collected as elements of a row vector
 */
Eigen::MatrixXd transformOneElectronIntegrals(const Eigen::MatrixXd& h, const Eigen::MatrixXd& T) {
    return T.adjoint() * h * T;
}


/** Given:
 *      - a rank-four tensor g, which contains the two-electron integrals in some basis B
 *      - a transformation matrix which represents the transformation of B into B'
 *
 * return the tensor g' in the new basis B'
 *
 *
 *
 *
 *  Note that the transformation matrix T is used as
 *
 *      B' = B T ,
 *
 *  where the basis vectors are collected as elements of a row vector
 */
Eigen::Tensor<double, 4> transformTwoElectronIntegrals(const Eigen::Tensor<double, 4>& g, const Eigen::MatrixXd& T) {

    // Since we're only getting T as a matrix, we should make the appropriate tensor to perform contractions
    // For the const Eigen::MatrixXd& argument, we need the const double in the template
    //      For more info, see: https://stackoverflow.com/questions/45283468/eigen-const-tensormap
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> T_tensor (T.data(), T.rows(), T.cols());


    // We will have to do four single contractions, so we specify the contraction indices
    // Eigen3 does not document its tensor contraction clearly, so see the accepted answer on stackoverflow (https://stackoverflow.com/a/47558349/7930415):
    //      Eigen3 does not accept a way to specify the output axes: instead, it retains the order from left to right of the axes that survive the contraction.
    //      This means that, in order to get the right ordering of the axes, we will have to swap axes

    // g(T U V W)  T^*(V R) -> a(T U R W) but we get a(T U W R)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair1 = {Eigen::IndexPair<int>(2, 0)};
    Eigen::array<int, 4> shuffle_1 {0, 1, 3, 2};

    // a(T U R W)  T(W S) -> b(T U R S) and we get b(T U R S), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair2 = {Eigen::IndexPair<int>(3, 0)};

    // T(U Q)  b(T U R S) -> c(T Q R S) but we get c(Q T R S)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair3 = {Eigen::IndexPair<int>(0, 1)};
    Eigen::array<int, 4> shuffle_3 {1, 0, 2, 3};

    // T^*(T P)  c(T Q R S) -> g'(P Q R S) and we get g_SO(p q r s), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair4 = {Eigen::IndexPair<int>(0, 0)};


    // Calculate the contractions. We write this as one large contraction to
    //  1) avoid storing intermediate contractions
    //  2) let Eigen3 figure out some optimizations
    Eigen::Tensor<double, 4> g_transformed = T_tensor.conjugate().contract(T_tensor.contract(g.contract(T_tensor.conjugate(), contraction_pair1).shuffle(shuffle_1).contract(T_tensor, contraction_pair2), contraction_pair3).shuffle(shuffle_3), contraction_pair4);

    return g_transformed;
};


/*
 *  AO AND SO CONVERSION WRAPPERS
 */

/** Given:
 *      - a matrix representation in an AO basis (f_AO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the SO basis
 */
Eigen::MatrixXd transform_AO_to_SO(const Eigen::MatrixXd& f_AO, const Eigen::MatrixXd& C) {
    return transformOneElectronIntegrals(f_AO, C);
}


/** Given:
 *      - a matrix representation in an SO basis (f_SO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the AO basis
 */
Eigen::MatrixXd transform_SO_to_AO(const Eigen::MatrixXd& f_SO, const Eigen::MatrixXd& C){
    Eigen::MatrixXd C_inverse = C.inverse();
    return transformOneElectronIntegrals(f_SO, C_inverse);
}


/** Given:
 *      - a rank-four tensor of two-electron integrals in an AO basis
 *      - an SO coefficient matrix (every column represents a spatial orbital)
 *
 *  transform and return the two-electron integrals in the SO basis
 */
Eigen::Tensor<double, 4> transform_AO_to_SO(const Eigen::Tensor<double, 4>& g_AO, const Eigen::MatrixXd& C) {
    return transformTwoElectronIntegrals(g_AO, C);
};


/*
 *  JACOBI ROTATIONS AND WRAPPERS
 */

/** Give the M-dimensional Jacobi rotation matrix (with an angle @param: theta in radians) for the orbitals P and Q (P < Q).
 *
 * M is the actual dimension of the matrix that is returned
 * P and Q represent the rows and columns, i.e. they start at 0
 *
 * Note that we work with the (cos, sin, -sin, cos) definition
 */
Eigen::MatrixXd jacobiRotationMatrix(size_t p, size_t q, double theta, size_t M) {

    if (p >= q) {
        throw std::invalid_argument("p should be smaller than q");
    }

    if ((M < p + 1) || (M < q + 1)) {
        throw std::invalid_argument("M should be larger than (p+1) and larger than (q+1).");
    }

    // The union of these two conditions also excludes M < 2


    double c = std::cos(theta);
    double s = std::sin(theta);

    // We'll start the construction with an identity matrix
    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(M, M);

    // And apply the Jacobi rotation as J = I * jacobi_rotation (cfr. B' = B T)
    J.applyOnTheRight(p, q, Eigen::JacobiRotation<double> (c, s));
    return J;
}


/**
 *  Using a Jacobi rotation with angle theta of the orbitals p and q, return the transformed one-electron integrals.
 *  This function is implemented using Eigen's Jacobi module.
 */
Eigen::MatrixXd rotateOneElectronIntegralsJacobi(const Eigen::MatrixXd& h, size_t p, size_t q, double theta) {

    double c = std::cos(theta);
    double s = std::sin(theta);

    // Initialize the rotated matrix by making a copy of the original matrix
    Eigen::MatrixXd h_rotated = h;

    // Use Eigen's Jacobi module to apply the Jacobi rotations directly (cfr. T.adjoint() * h * T)
    Eigen::JacobiRotation<double> jacobi (c, s);
    h_rotated.applyOnTheLeft(p, q, jacobi.adjoint());
    h_rotated.applyOnTheRight(p, q, jacobi);

    return h_rotated;
}


/**
 *  Using a Jacobi rotation with angle theta of the orbitals p and q, return the transformed two-electron integrals.
 *  While waiting for an analogous Eigen::Tensor Jacobi module, this function is just a wrapper around transformTwoElectronIntegrals using a Jacobi rotation matrix.
 */
Eigen::Tensor<double, 4> rotateTwoElectronIntegralsJacobi(const Eigen::Tensor<double, 4>& g, size_t p, size_t q, double theta) {

    auto dim = static_cast<size_t>(g.dimension(1));  // g.dimension() returns a long
    Eigen::MatrixXd J = jacobiRotationMatrix(p, q, theta, dim);

    return transformTwoElectronIntegrals(g, J);
};


}  // namespace transformations
}  // namespace libwint
