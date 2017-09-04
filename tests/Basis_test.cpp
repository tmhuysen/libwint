#define BOOST_ALL_DYN_LINK

#define BOOST_TEST_MODULE "Molecule class"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain
#include "Molecule.hpp"
#include "Basis.hpp"
#include <libint2.hpp>


/* Reads a
 *
 */
void read_array_from_file(const std::string &filename, Eigen::MatrixXf& M){
    std::ifstream file (filename);

    if (file.is_open()) {
        std::string line;
        while (std::getline (file, line)) {
            std::istringstream is (line);

            int i;
            int j;
            float value;

            is >> i >> j >> value;
            M(i, j) = value;
        }

        file.close();
    }
}


void read_array_from_file(const std::string &filename, Eigen::Tensor<double, 4>& M){
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


void check_equal_arrays(const Eigen::MatrixXf& M, const std::string& filename){
    auto dim = M.rows();    // The given matrices are symmetric so M.rows() == M.cols() is the dimension of the matrix
                            // This is also the dimension of T, V and the rank-4 tensor tei
    Eigen::MatrixXf M_data (dim, dim);

    read_array_from_file(filename, M_data);
    BOOST_CHECK(M.isApprox(M_data));
}

void check_equal_arrays(Eigen::Tensor<double, 4>& M, const std::string& filename){
    auto dim = M.NumIndices;    // The given matrices are symmetric so M.rows() == M.cols() is the dimension of the matrix
                                // This is also the dimension of T, V and the rank-4 tensor tei
    Eigen::Tensor<double, 4> M_data (dim, dim, dim, dim);

    auto tolerance = 0.00000001;    // Hard-coded tolerance of 1e-08

    // Since Eigen::Tensor doesn't have an isApprox yet, we will check every pair of values manually
    read_array_from_file(filename, M_data);
    for (int i = 0; i < dim ; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    // FIXME: add documentation for chemist vs physicist integrals

                    BOOST_CHECK_LT(M_data(i,j,k,l) - M(i,j,k,l), tolerance);
                }
            }
        }
    }
}



BOOST_AUTO_TEST_CASE( constructor ){
    // Initialize libint2
    libint2::initialize();

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    std::string basis_name = "STO-3G";
    Molecule water (xyzfilename);

    Basis basis (water, basis_name);


    BOOST_CHECK_EQUAL(basis.molecule.natoms(), 3);
    BOOST_CHECK_EQUAL(basis.name, "STO-3G");

    // Finalize libint2
    libint2::finalize();
}

BOOST_AUTO_TEST_CASE( integrals ){
    // Initialize libint2
    libint2::initialize();

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    std::string basis_name = "STO-3G";
    Molecule water (xyzfilename);

    Basis basis (water, basis_name);

    auto S = basis.compute_overlap_integrals();
    auto T = basis.compute_kinetic_integrals();
    auto V = basis.compute_nuclear_integrals();

    auto tei = basis.compute_two_electron_integrals();


    const auto overlap_data = "/Users/laurentlemmens/Software/libint-eigen/tests/ref_data/overlap.data";
    const auto kinetic_data = "/Users/laurentlemmens/Software/libint-eigen/tests/ref_data/kinetic.data";
    const auto nuclear_data = "/Users/laurentlemmens/Software/libint-eigen/tests/ref_data/nuclear.data";
    const auto two_electron_data = "/Users/laurentlemmens/Software/libint-eigen/tests/ref_data/two_electron.data";

    check_equal_arrays(S, overlap_data);
    check_equal_arrays(T, kinetic_data);
    check_equal_arrays(V, nuclear_data);
    check_equal_arrays(tei, two_electron_data);


    // Finalize libint2
    libint2::finalize();
}


