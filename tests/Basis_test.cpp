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


BOOST_AUTO_TEST_CASE( constructor ) {
    // Initialize libint2
    libint2::initialize();

    const auto xyzfilename = "../../docs/h2o.xyz";
    const auto basis_name = "STO-3G";

    Wrapper::Molecule water (xyzfilename);
    Wrapper::Basis basis (water, basis_name);

    BOOST_CHECK_EQUAL(basis.name, "STO-3G");
    BOOST_CHECK_EQUAL(basis.nbf(), 7);

    // Finalize libint2
    libint2::finalize();
}
/*

BOOST_AUTO_TEST_CASE( horton_integrals_h2o_sto3g ) {
    // Initialize libint2
    libint2::initialize();

    std::string xyzfilename = "../../docs/h2o.xyz";
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


BOOST_AUTO_TEST_CASE( szabo_h2_sto3g ) {
    // We will follow the example in Szabo, section 3.5.2, where it is stated that R = 1.4 a.u. = 0.740848 Angstrom
    libint2::initialize();

    // Specify the data
    std::string xyzfilename = "../../docs/h2.xyz";
    std::string basis_name = "STO-3G";

    // Create a Molecule and a Basis
    Wrapper::Molecule h2 (xyzfilename);
    BOOST_CHECK(std::abs(h2.atoms[1].z - 1.4) < 1.0e-6);    // Check if the conversion from Angstrom to a.u. is correct
    Wrapper::Basis basis (h2, basis_name);
    BOOST_CHECK_EQUAL(basis.nbf(), 2);                      // Check if there are only two basis functions

    // Calculate S, T, V and H_core
    Eigen::MatrixXd S = basis.compute_overlap_integrals();
    Eigen::MatrixXd T = basis.compute_kinetic_integrals();
    Eigen::MatrixXd V = basis.compute_nuclear_integrals();
    Eigen::MatrixXd H_core = T + V;

    // Fill in the reference values from Szabo
    Eigen::MatrixXd S_ref (2, 2);
    S_ref << 1.0, 0.6593,
             0.6593, 1.0;

    Eigen::MatrixXd T_ref (2, 2);
    T_ref << 0.7600, 0.2365,
             0.2365, 0.7600;

    Eigen::MatrixXd H_core_ref (2, 2);
    H_core_ref << -1.1204, -0.9584,
                  -0.9584, -1.1204;

    BOOST_CHECK(S.isApprox(S_ref, 1.0e-4));
    BOOST_CHECK(T.isApprox(T_ref, 1.0e-4));
    BOOST_CHECK(H_core.isApprox(H_core_ref, 1.0e-4));


    // Calculate the two-electron integrals, and check the unique values listed in Szabo. These are given in chemist's notation in Szabo, so this confirms that this wrapper gives them in chemist's notation as well.
    Eigen::Tensor<double, 4> tei = basis.compute_two_electron_integrals();

    BOOST_CHECK(std::abs(tei(0,0,0,0) - 0.7746) < 1.0e-4);
    BOOST_CHECK(std::abs(tei(0,0,0,0) - tei(1,1,1,1)) < 1.0e-12);

    BOOST_CHECK(std::abs(tei(0,0,1,1) - 0.5697) < 1.0e-4);

    BOOST_CHECK(std::abs(tei(1,0,0,0) - 0.4441) < 1.0e-4);
    BOOST_CHECK(std::abs(tei(1,0,0,0) - tei(1,1,1,0)) < 1.0e-12);

    BOOST_CHECK(std::abs(tei(1,0,1,0) - 0.2970) < 1.0e-4);


    libint2::finalize();
}