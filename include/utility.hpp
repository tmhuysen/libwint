#ifndef LIBWINT_UTILITY_HPP
#define LIBWINT_UTILITY_HPP

#include <Eigen/Dense>
#include <string>
#include <unsupported/Eigen/CXX11/Tensor>


namespace libwint {
namespace utility {

/** Read an array from a given filename line by line, and add the elements to a given matrix (rank-2 tensor)
*/
void read_array_from_file(const std::string& filename, Eigen::MatrixXd& M);


/** Read an array from a given filename line by line, and add the elements to a given tensor (rank-4 tensor)
*/
void read_array_from_file(const std::string& filename, Eigen::Tensor<double, 4>& M);


/** Return if two rank-4 tensors are approximately equal
 */
bool are_equal(const Eigen::Tensor<double, 4>& M, const Eigen::Tensor<double, 4>& T, const double tolerance);


/** Print the contents of a rank-four tensor in a fashionable way
 */
void print(const Eigen::Tensor<double, 4>& T);

}  // namespace utility
}  // namespace libwint


#endif // LIBWINT_UTILITY_HPP
