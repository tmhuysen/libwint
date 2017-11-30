#ifndef LIBWRP_TRANSFORMATIONS_HPP
#define LIBWRP_TRANSFORMATIONS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


namespace libwrp {


/** Given:
 *      - a matrix of one-electron integrals in an AO basis
 *      - an SO coefficient matrix (every column represents a spatial orbital)
 *
 *  transform and return the one-electron integrals in the SO basis
 */
Eigen::MatrixXd transform_AO_to_SO(Eigen::MatrixXd &one_electron_integrals, Eigen::MatrixXd &C);


/** Given:
 *      - a rank-four tensor of two-electron integrals in an AO basis
 *      - an SO coefficient matrix (every column represents a spatial orbital)
 *
 *  transform and return the two-electron integrals in the SO basis
 */
Eigen::Tensor<double, 4> transform_AO_to_SO(Eigen::Tensor<double, 4> &g_AO, Eigen::MatrixXd &C);


}  // namespace libwrp

#endif // LIBWRP_TRANSFORMATIONS_HPP
