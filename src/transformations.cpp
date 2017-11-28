#include "transformations.hpp"

#include <iostream>

/** Given:
 *      - a matrix of one-electron integrals in an AO basis
 *      - an SO coefficient matrix (every column represents a spatial orbital)
 *
 *  transform and return the one-electron integrals in the SO basis
 */
Eigen::MatrixXd libwrp::transform_AO_integrals_to_SO(Eigen::MatrixXd& one_electron_integrals, Eigen::MatrixXd& C) {
    return C.adjoint() * one_electron_integrals * C;
}


/** Given:
 *      - a rank-four tensor of two-electron integrals in an AO basis
 *      - an SO coefficient matrix (every column represents a spatial orbital)
 *
 *  transform and return the two-electron integrals in the SO basis
 */
Eigen::Tensor<double, 4> libwrp::transform_AO_integrals_to_SO(Eigen::Tensor<double, 4>& two_electron_integrals, Eigen::MatrixXd& C) {

//    // Since we're only getting C as a matrix, we should make the appropriate tensor to perform contractions
//    Eigen::TensorMap<Eigen::Tensor<double, 2>> C_tensor (C.data(), C.rows(), C.cols());
//
//
//    // We will have to do four single contractions. Specify the contraction indices:
//    //  g(mu nu lambda rho)  C^*(lambda r)
//    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair1 = {Eigen::IndexPair<int>(2, 0)};
//
//    //  a(mu nu r rho)  C(rho s)
//    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair2 = {Eigen::IndexPair<int>(3, 0)};
//
//    //  C(nu q)  b(mu nu r s)
//    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair3 = {Eigen::IndexPair<int>(0, 1)};
//
//    //  C^*(mu p)  c(mu q r s)
//    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair4 = {Eigen::IndexPair<int>(0, 0)};
//
//
//    // Calculate the contractions
//    Eigen::Tensor<double, 4> a = two_electron_integrals.contract(C_tensor.conjugate(), contraction_pair1);
//    Eigen::Tensor<double, 4> b = a.contract(C_tensor, contraction_pair2);
//    Eigen::Tensor<double, 4> c = C_tensor.contract(b, contraction_pair3);
//    Eigen::Tensor<double, 4> g_SO = C_tensor.conjugate().contract(b, contraction_pair4);
//
//    return g_SO;

    Eigen::MatrixXd C_conjugate = C.conjugate();  // Compute this first because we don't want to recompute this in the loops


    auto dim = static_cast<size_t>(C.cols());  // Eigen3 uses an int for cols
    Eigen::Tensor<double, 4> g_SO (dim, dim, dim, dim);
    g_SO.setZero();

    for (size_t p = 0; p < dim; p++) {
        for (size_t q = 0; q < dim; q++) {
            for (size_t r = 0; r < dim; r++) {
                for (size_t s = 0; s < dim; s++) {

                    for (size_t mu = 0; mu < dim; mu++) {
                        for (size_t nu = 0; nu < dim; nu++) {
                            for (size_t lambda = 0; lambda < dim; lambda++) {
                                for (size_t rho = 0; rho < dim; rho++) {

                                    g_SO(p, q, r, s) += C_conjugate(mu, p) * C(nu, q) * two_electron_integrals(mu, nu, lambda, rho) * C_conjugate(lambda, r) * C(rho, s);

                                }
                            }
                        }
                    }  // mu nu lambda rho

                }
            }
        }
    }  // pqrs

    std::cout << g_SO << std::endl;


    return g_SO;
};
