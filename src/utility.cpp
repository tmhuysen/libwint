#include "utility.hpp"


/** Read an array from a given filename line by line, and add the elements to a given matrix (rank-2 tensor)
*/
void libwrp::utility::read_array_from_file(const std::string& filename, Eigen::MatrixXd& M){
    std::ifstream file (filename);

    if (file.is_open()) {
        std::string line;
        while (std::getline (file, line)) {
            std::istringstream is (line);

            int i;
            int j;
            double value;

            is >> i >> j >> value;
            M(i, j) = value;
        }

        file.close();
    }
}


/** Read an array from a given filename line by line, and add the elements to a given tensor (rank-4 tensor)
*/
void libwrp::utility::read_array_from_file(const std::string& filename, Eigen::Tensor<double, 4>& M){
    std::ifstream file (filename);

    if (file.is_open()) {
        std::string line;
        while (std::getline (file, line)) {
            std::istringstream is (line);

            int i;
            int j;
            int k;
            int l;
            float value;

            is >> i >> j >> k >> l >> value;
            M(i, j, k, l) = value;
            std::cout << i << ' ' << j << ' ' << k << ' ' << l << "  " << M(i, j, k, l) << std::endl;
        }

        file.close();
    }
}


/** Return if two rank-4 tensors are approximately equal
 */
bool libwrp::utility::are_equal(const Eigen::Tensor<double, 4>& M, const Eigen::Tensor<double, 4>& T, const double tolerance){
    auto dim = M.dimension(0);

    // Since Eigen::Tensor doesn't have an isApprox yet, we will check every pair of values manually
    for (int i = 0; i < dim ; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    if (std::abs(M(i,j,k,l) - T(i,j,k,l)) > tolerance) {
                        std::cout << M(i,j,k,l) << ' ' << T(i,j,k,l) << std::endl;
                        return false;
                    }
                }
            }
        }
    } // rank-4 tensor traversing
    return true;
}


/** Print the contents of a rank-four tensor in a fashionable way
 */
void libwrp::utility::print(const Eigen::Tensor<double, 4>& T) {
    auto dim = T.dimension(0);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    std::cout << i << ' ' << j << ' ' << k << ' ' << l << "  " << T(i, j, k, l) << std::endl;
                }
            }
        }
    }
}
