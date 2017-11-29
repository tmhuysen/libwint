#define BOOST_TEST_MODULE "Basis"

#include "Basis.hpp"
#include "utility.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE( constructor ) {
    // Initialize libint2
    libint2::initialize();

    const std::string xyzfilename = "../tests/ref_data/h2o.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    const std::string basis_name = "STO-3G";

    libwrp::Molecule water (xyzfilename);
    libwrp::Basis basis (water, basis_name);

    BOOST_CHECK_EQUAL(basis.name, "STO-3G");
    BOOST_CHECK_EQUAL(basis.nbf(), 7);

    // Test if a reference to the Molecule object is actually made
    BOOST_CHECK_EQUAL(&basis.molecule, &water);

    // Finalize libint2
    libint2::finalize();
}


BOOST_AUTO_TEST_CASE( horton_integrals_h2o_sto3g ) {
    // Initialize libint2
    libint2::initialize();

    const std::string xyzfilename = "../tests/ref_data/h2o.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    const std::string basis_name = "STO-3G";
    libwrp::Molecule water (xyzfilename);
    libwrp::Basis basis (water, basis_name);
    auto nbf = basis.nbf();

    basis.compute_integrals();

    Eigen::MatrixXd S_test (nbf, nbf);
    Eigen::MatrixXd T_test (nbf, nbf);
    Eigen::MatrixXd V_test (nbf, nbf);
    Eigen::Tensor<double, 4> tei_test (nbf, nbf, nbf, nbf);

    libwrp::utility::read_array_from_file("../tests/ref_data/overlap.data", S_test);
    libwrp::utility::read_array_from_file("../tests/ref_data/kinetic.data", T_test);
    libwrp::utility::read_array_from_file("../tests/ref_data/nuclear.data", V_test);
    libwrp::utility::read_array_from_file("../tests/ref_data/two_electron.data", tei_test);

    BOOST_CHECK(basis.S.isApprox(S_test, 1.0e-8));
    BOOST_CHECK(basis.T.isApprox(T_test, 1.0e-8));
    BOOST_CHECK(basis.V.isApprox(V_test, 1.0e-8));
    BOOST_CHECK(libwrp::utility::are_equal(basis.tei, tei_test, 1.0e-6));

    // Finalize libint2
    libint2::finalize();
}


BOOST_AUTO_TEST_CASE( szabo_h2_sto3g ) {
    // We will follow the example in Szabo, section 3.5.2, where it is stated that R = 1.4 a.u. = 0.740848 Angstrom
    libint2::initialize();

    // Specify the data
    const std::string xyzfilename = "../tests/ref_data/h2.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    const std::string basis_name = "STO-3G";

    // Create a Molecule and a Basis
    libwrp::Molecule h2 (xyzfilename);
    BOOST_CHECK(std::abs(h2.atoms[1].z - 1.4) < 1.0e-6);    // Check if the conversion from Angstrom to a.u. is correct
    libwrp::Basis basis (h2, basis_name);
    BOOST_CHECK_EQUAL(basis.nbf(), 2);                      // Check if there are only two basis functions

    // Calculate S, T, V and H_core
    basis.compute_integrals();
    Eigen::MatrixXd H_core = basis.T + basis.V;

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

    BOOST_CHECK(basis.S.isApprox(S_ref, 1.0e-4));
    BOOST_CHECK(basis.T.isApprox(T_ref, 1.0e-4));
    BOOST_CHECK(H_core.isApprox(H_core_ref, 1.0e-4));


    // Calculate the two-electron integrals, and check the unique values listed in Szabo. These are given in chemist's notation in Szabo, so this confirms that this libwrp gives them in chemist's notation as well.
    BOOST_CHECK(std::abs(basis.tei(0,0,0,0) - 0.7746) < 1.0e-4);
    BOOST_CHECK(std::abs(basis.tei(0,0,0,0) - basis.tei(1,1,1,1)) < 1.0e-12);

    BOOST_CHECK(std::abs(basis.tei(0,0,1,1) - 0.5697) < 1.0e-4);

    BOOST_CHECK(std::abs(basis.tei(1,0,0,0) - 0.4441) < 1.0e-4);
    BOOST_CHECK(std::abs(basis.tei(1,0,0,0) - basis.tei(1,1,1,0)) < 1.0e-12);

    BOOST_CHECK(std::abs(basis.tei(1,0,1,0) - 0.2970) < 1.0e-4);

    libint2::finalize();
}
