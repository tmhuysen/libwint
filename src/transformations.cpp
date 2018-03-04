#include "transformations.hpp"

#include <iostream>


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
Eigen::MatrixXd libwint::transformOneElectronIntegrals(Eigen::MatrixXd& h, Eigen::MatrixXd& T) {
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
Eigen::Tensor<double, 4> libwint::transform_two_electron_integrals(Eigen::Tensor<double, 4>& g, Eigen::MatrixXd& T) {

    // Since we're only getting T as a matrix, we should make the appropriate tensor to perform contractions
    Eigen::TensorMap<Eigen::Tensor<double, 2>> T_tensor (T.data(), T.rows(), T.cols());


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



/** Given:
 *      - a matrix representation in an AO basis (f_AO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the SO basis
 */
Eigen::MatrixXd libwint::transform_AO_to_SO(Eigen::MatrixXd& f_AO, Eigen::MatrixXd& C) {
    return transformOneElectronIntegrals(f_AO, C);
}


/** Given:
 *      - a matrix representation in an SO basis (f_SO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the AO basis
 */
Eigen::MatrixXd libwint::transform_SO_to_AO(Eigen::MatrixXd& f_SO, Eigen::MatrixXd& C){
    Eigen::MatrixXd C_inverse = C.inverse();
    return transformOneElectronIntegrals(f_SO, C_inverse);
}


/** Given:
 *      - a rank-four tensor of two-electron integrals in an AO basis
 *      - an SO coefficient matrix (every column represents a spatial orbital)
 *
 *  transform and return the two-electron integrals in the SO basis
 */
Eigen::Tensor<double, 4> libwint::transform_AO_to_SO(Eigen::Tensor<double, 4>& g_AO, Eigen::MatrixXd& C) {
    return transform_two_electron_integrals(g_AO, C);
};


/** Given a unitary matrix U that transforms a basis B into B', return the one-electron integrals in the rotated basis
 *
 * Note that the basis transformation is explicitly written as (B' = B U)
 */
Eigen::MatrixXd libwint::rotate_integrals(Eigen::MatrixXd& h, Eigen::MatrixXd& U) {

    // Check if the given matrix U is unitary
    if (!U.isUnitary()) {
        throw std::invalid_argument("The given matrix U is not unitary.");
    }

    return transformOneElectronIntegrals(h, U);
}


/** Given a unitary matrix U that transforms a basis B into B', return the two-electron integrals in the rotated basis
 *
 * Note that the basis transformation is explicitly written as (B' = B U)
 */
Eigen::Tensor<double, 4> libwint::rotate_integrals(Eigen::Tensor<double, 4>& g, Eigen::MatrixXd& U) {

    // Check if the given matrix U is unitary
    if (!U.isUnitary()) {
        throw std::invalid_argument("The given matrix U is not unitary.");
    }

    return transform_two_electron_integrals(g, U);
};


/** Give the M-dimensional Jacobi rotation matrix (with an angle theta) for the orbitals P and Q (P < Q).
 *
 * M is the actual dimension of the matrix that is returned
 * P and Q represent the rows and columns, i.e. they start at 0
 *
 * Note that we work with the (cos, sin, -sin, cos) definition
 */
Eigen::MatrixXd libwint::jacobi_rotation_matrix(size_t P, size_t Q, double theta, size_t M) {

    if (P >= Q) {
        throw std::invalid_argument("P should be smaller than Q");
    }

    if ((M < P + 1) || (M < Q + 1)) {
        throw std::invalid_argument("M should be larger than (P+1) and larger than (Q+1).");
    }

    // The union of these two conditions also excludes M < 2



    // We'll start the construction with an identity matrix
    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(M, M);

    // Add the Jacobi rotation terms
    J(P, P) = std::cos(theta);
    J(P, Q) = std::sin(theta);
    J(Q, P) = -std::sin(theta);
    J(Q, Q) = std::cos(theta);

    return J;
}


/** Using a Jacobi rotation with angle theta of the orbitals P and Q, return the transformed one-electron integrals.
 *
 * In this function, I've implemented it so it can be checked to be correct.
 * In the analytical derivation, I have explicitly assumed that we are working with a symmetric matrix h (h_PQ = h_QP)
 */
Eigen::MatrixXd libwint::rotate_one_electron_integrals_jacobi(Eigen::MatrixXd& h, size_t P, size_t Q, double theta) {

    // Assert the assumption of a symmetric matrix
    assert(h.isApprox(h.transpose()));

    // Initialize the rotated matrix by making a copy of the original matrix
    Eigen::MatrixXd h_rotated = h;

    // Since we have a Jacobi rotation, we can directly fill in rows and columns P and Q
    // Update the P-th and Q-th row
    for (size_t S = 0; S < h.cols(); S++) {
        h_rotated(P,S) += h(P,S) * (std::cos(theta) - 1) - h(Q,S) * std::sin(theta);
        h_rotated(Q,S) += h(Q,S) * (std::cos(theta) - 1) + h(P,S) * std::sin(theta);
    }

    // Update the P-th and Q-th column
    for (size_t R = 0; R < h.rows(); R++) {
        h_rotated(R,P) += h(R,P) * (std::cos(theta) - 1) - h(R,Q) * std::sin(theta);
        h_rotated(R,Q) += h(R,Q) * (std::cos(theta) - 1) + h(R,P) * std::sin(theta);
    }

    // Update the four intersections (P,P) (P,Q) (Q,P) (Q,Q)
    h_rotated(P,P) += h(P,P) * std::pow(std::cos(theta) - 1, 2) - 2 * h(P,Q) * (std::cos(theta) - 1) * std::sin(theta) + h(Q,Q) * std::pow(std::sin(theta), 2);
    h_rotated(P,Q) += h(P,Q) * std::pow(std::cos(theta) - 1, 2) + (h(P,P) - h(Q,Q)) * (std::cos(theta) - 1) * std::sin(theta) - h(P,Q) * std::pow(std::sin(theta), 2);
    h_rotated(Q,P) += h(P,Q) * std::pow(std::cos(theta) - 1, 2) + (h(P,P) - h(Q,Q)) * (std::cos(theta) - 1) * std::sin(theta) - h(P,Q) * std::pow(std::sin(theta), 2);
    h_rotated(Q,Q) += h(Q,Q) * std::pow(std::cos(theta) - 1, 2) + 2 * h(P,Q) * (std::cos(theta) - 1) * std::sin(theta) + h(P,P) * std::pow(std::sin(theta), 2);

    return h_rotated;
}
