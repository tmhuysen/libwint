#define BOOST_ALL_DYN_LINK

#define BOOST_TEST_MODULE "Basis"

#include "Basis.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

/** Read an array from a given filename line by line, and add the elements to a given matrix (rank-2 tensor)
*/
void read_array_from_file(const std::string& filename, Eigen::MatrixXd& M){
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
void read_array_from_file(const std::string& filename, Eigen::Tensor<double, 4>& M){
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
        }

        file.close();
    }
}

/** Return if two rank-4 tensors are approximately equal
 */
bool are_equal(const Eigen::Tensor<double, 4>& M, const Eigen::Tensor<double, 4>& T, const double tolerance){
    auto dim = M.NumIndices;

    // Since Eigen::Tensor doesn't have an isApprox yet, we will check every pair of values manually
    for (int i = 0; i < dim ; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    if (std::abs(M(i,j,k,l) - T(i,j,k,l)) > tolerance) {
                        return false;
                    }
                }
            }
        }
    } // rank-4 tensor traversing
    return true;
}


BOOST_AUTO_TEST_CASE( constructor ){
    // Initialize libint2
    libint2::initialize();

    const auto xyzfilename = "../../docs/h2o.xyz";
    std::string basis_name = "STO-3G";
    Wrapper::Molecule water (xyzfilename);
    Wrapper::Basis basis (water, basis_name);

    BOOST_CHECK_EQUAL(basis.name, "STO-3G");
    BOOST_CHECK_EQUAL(basis.nbf(), 7);

    // Finalize libint2
    libint2::finalize();
}


BOOST_AUTO_TEST_CASE( integrals ){
    // Initialize libint2
    libint2::initialize();

    const auto xyzfilename = "../../docs/h2o.xyz";
    std::string basis_name = "STO-3G";
    Wrapper::Molecule water (xyzfilename);
    Wrapper::Basis basis (water, basis_name);
    auto nbf = basis.nbf();

    Eigen::MatrixXd S = basis.compute_overlap_integrals();
    Eigen::MatrixXd T = basis.compute_kinetic_integrals();
    Eigen::MatrixXd V = basis.compute_nuclear_integrals();
    Eigen::Tensor<double, 4> tei = basis.compute_two_electron_integrals();

    Eigen::MatrixXd S_test (nbf, nbf);
    Eigen::MatrixXd T_test (nbf, nbf);
    Eigen::MatrixXd V_test (nbf, nbf);
    Eigen::Tensor<double, 4> tei_test (nbf, nbf, nbf, nbf);

    read_array_from_file("../../tests/ref_data/overlap.data", S_test);
    read_array_from_file("../../tests/ref_data/kinetic.data", T_test);
    read_array_from_file("../../tests/ref_data/nuclear.data", V_test);
    read_array_from_file("../../tests/ref_data/two_electron.data", tei_test);

    BOOST_CHECK(S.isApprox(S_test, 1.0e-8));
    BOOST_CHECK(T.isApprox(T_test, 1.0e-8));
    BOOST_CHECK(V.isApprox(V_test, 1.0e-8));
    BOOST_CHECK(are_equal(tei, tei_test, 1.0e-6));

    // Finalize libint2
    libint2::finalize();
}
