#include "transformations.hpp"


/** Given:
 *      - a matrix representation in an AO basis (f_AO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the SO basis
 */
Eigen::MatrixXd libwrp::transform_AO_to_SO(Eigen::MatrixXd& f_AO, Eigen::MatrixXd& C) {
    return C.adjoint() * f_AO * C;
}


/** Given:
 *      - a matrix representation in an SO basis (f_SO)
 *      - an SO coefficient matrix (every column represents a spatial orbital) (C)
 *
 *  transform and return the matrix in the AO basis
 */
Eigen::MatrixXd libwrp::transform_SO_to_AO(Eigen::MatrixXd& f_SO, Eigen::MatrixXd& C){
    return C.inverse().adjoint() * f_SO * C.inverse();
}


/** Given:
 *      - a rank-four tensor of two-electron integrals in an AO basis
 *      - an SO coefficient matrix (every column represents a spatial orbital)
 *
 *  transform and return the two-electron integrals in the SO basis
 */
Eigen::Tensor<double, 4> libwrp::transform_AO_to_SO(Eigen::Tensor<double, 4>& g_AO, Eigen::MatrixXd& C) {

    // Since we're only getting C as a matrix, we should make the appropriate tensor to perform contractions
    Eigen::TensorMap<Eigen::Tensor<double, 2>> C_tensor (C.data(), C.rows(), C.cols());


    // We will have to do four single contractions, so we specify the contraction indices
    // Eigen3 does not document its tensor contraction clearly, so see the accepted answer on stackoverflow (https://stackoverflow.com/a/47558349/7930415):
    //      Eigen3 does not accept a way to specify the output axes: instead, it retains the order from left to right of the axes that survive the contraction.
    //      This means that, in order to get the right ordering of the axes, we will have to swap axes

    // g_AO(mu nu lambda rho)  C^*(lambda r) -> a(mu nu r rho) but we get a(mu nu rho lambda)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair1 = {Eigen::IndexPair<int>(2, 0)};
    Eigen::array<int, 4> shuffle_1 {0, 1, 3, 2};

    // a(mu nu r rho)  C(rho s) -> b(mu nu r s) and we get b(mu nu r s), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair2 = {Eigen::IndexPair<int>(3, 0)};

    // C(nu q)  b(mu nu r s) -> c(mu q r s) but we get c(q mu r s)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair3 = {Eigen::IndexPair<int>(0, 1)};
    Eigen::array<int, 4> shuffle_3 {1, 0, 2, 3};

    // C^*(mu p)  c(mu q r s) -> g_SO(p q r s) and we get g_SO(p q r s), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair4 = {Eigen::IndexPair<int>(0, 0)};


    // Calculate the contractions. We write this as one large contraction to
    //  1) avoid storing intermediate contractions
    //  2) let Eigen3 figure out some optimizations
    Eigen::Tensor<double, 4> g_SO = C_tensor.conjugate().contract(C_tensor.contract(g_AO.contract(C_tensor.conjugate(), contraction_pair1).shuffle(shuffle_1).contract(C_tensor, contraction_pair2), contraction_pair3).shuffle(shuffle_3), contraction_pair4);

    return g_SO;
};
