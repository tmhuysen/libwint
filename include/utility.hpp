#ifndef LIBWRP_UTILITY_HPP
#define LIBWRP_UTILITY_HPP


namespace libwrp::utility {

/** Read an array from a given filename line by line, and add the elements to a given matrix (rank-2 tensor)
*/
void libwrp::utility::read_array_from_file(const std::string& filename, Eigen::MatrixXd& M);


/** Read an array from a given filename line by line, and add the elements to a given tensor (rank-4 tensor)
*/
void libwrp::utility::read_array_from_file(const std::string& filename, Eigen::Tensor<double, 4>& M);


/** Return if two rank-4 tensors are approximately equal
 */
bool libwrp::utility::are_equal(const Eigen::Tensor<double, 4>& M, const Eigen::Tensor<double, 4>& T, const double tolerance);


/** Print the contents of a rank-four tensor in a fashionable way
 */
void libwrp::utility::print(const Eigen::Tensor<double, 4>& T);

}  // namespace libwrp::utility


#endif // LIBWRP_UTILITY_HPP
