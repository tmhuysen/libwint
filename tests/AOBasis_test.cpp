#define BOOST_TEST_MODULE "Basis"


#include "AOBasis.hpp"

#include "cpputil.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE( basis_constructor ) {

    // Check the number of basis functions in water
    libwint::Molecule water ("../tests/ref_data/h2o.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::AOBasis basis (water, "STO-3G");
    basis.calculateIntegrals();

    BOOST_CHECK_EQUAL(basis.calculateNumberOfBasisFunctions(), 7);
}


BOOST_AUTO_TEST_CASE( horton_integrals_h2o_sto3g ) {

    libwint::Molecule water ("../tests/ref_data/h2o.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::AOBasis basis (water, "STO-3G");
    basis.calculateIntegrals();
    size_t nbf = basis.calculateNumberOfBasisFunctions();

    // Read in reference data from Horton
    Eigen::MatrixXd ref_S (nbf, nbf);
    Eigen::MatrixXd ref_T (nbf, nbf);
    Eigen::MatrixXd ref_V (nbf, nbf);
    Eigen::Tensor<double, 4> ref_teri (nbf, nbf, nbf, nbf);

    cpputil::io::readArrayFromFile("../tests/ref_data/h2o_sto-3g_overlap.data", ref_S);
    cpputil::io::readArrayFromFile("../tests/ref_data/h2o_sto-3g_kinetic.data", ref_T);
    cpputil::io::readArrayFromFile("../tests/ref_data/h2o_sto-3g_nuclear.data", ref_V);
    cpputil::io::readArrayFromFile("../tests/ref_data/h2o_sto-3g_two_electron.data", ref_teri);


    BOOST_CHECK(basis.get_S().isApprox(ref_S, 1.0e-8));
    BOOST_CHECK(basis.get_T().isApprox(ref_T, 1.0e-8));
    BOOST_CHECK(basis.get_V().isApprox(ref_V, 1.0e-8));
    BOOST_CHECK(cpputil::linalg::areEqual(basis.get_g(), ref_teri, 1.0e-6));
}


BOOST_AUTO_TEST_CASE( szabo_h2_sto3g ) {

    // We will follow the example in Szabo, section 3.5.2, where it is stated that R = 1.4 a.u. = 0.740848 Angstrom
    libwint::Molecule h2 ("../tests/ref_data/h2.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::AOBasis basis (h2, "STO-3G");
    basis.calculateIntegrals();
    BOOST_CHECK_EQUAL(basis.calculateNumberOfBasisFunctions(), 2);

    // Calculate S, T, V and H_core
    Eigen::MatrixXd H_core = basis.get_T() + basis.get_V();

    // Fill in the reference values from Szabo
    Eigen::MatrixXd ref_S (2, 2);
    ref_S << 1.0, 0.6593,
             0.6593, 1.0;

    Eigen::MatrixXd ref_T (2, 2);
    ref_T << 0.7600, 0.2365,
             0.2365, 0.7600;

    Eigen::MatrixXd ref_H_core (2, 2);
    ref_H_core << -1.1204, -0.9584,
                  -0.9584, -1.1204;

    BOOST_CHECK(basis.get_S().isApprox(ref_S, 1.0e-4));
    BOOST_CHECK(basis.get_T().isApprox(ref_T, 1.0e-4));
    BOOST_CHECK(H_core.isApprox(ref_H_core, 1.0e-4));


    // Calculate the two-electron integrals, and check the unique values listed in Szabo.
    // These are given in chemist's notation in Szabo, so this confirms that this libwint gives them in chemist's notation as well.
    BOOST_CHECK(std::abs(basis.get_g()(0,0,0,0) - 0.7746) < 1.0e-4);
    BOOST_CHECK(std::abs(basis.get_g()(0,0,0,0) - basis.get_g()(1,1,1,1)) < 1.0e-12);

    BOOST_CHECK(std::abs(basis.get_g()(0,0,1,1) - 0.5697) < 1.0e-4);

    BOOST_CHECK(std::abs(basis.get_g()(1,0,0,0) - 0.4441) < 1.0e-4);
    BOOST_CHECK(std::abs(basis.get_g()(1,0,0,0) - basis.get_g()(1,1,1,0)) < 1.0e-12);

    BOOST_CHECK(std::abs(basis.get_g()(1,0,1,0) - 0.2970) < 1.0e-4);
}
