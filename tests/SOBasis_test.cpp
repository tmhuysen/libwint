#define BOOST_TEST_MODULE "SOBasis"


#include "SOBasis.hpp"

#include "transformations.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain



BOOST_AUTO_TEST_CASE ( transform_jacobi ) {

    // Create an SOBasis instance with a coefficient matrix being the identity matrix (little hack that we can use to test transformations)
    libwint::Molecule water ("../tests/ref_data/h2o.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::AOBasis ao_basis (water, "STO-3G");
    size_t K = ao_basis.calculateNumberOfBasisFunctions();
    ao_basis.calculateIntegrals();

    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(K, K);
    libwint::SOBasis so_basis (ao_basis, C);


    // Test if the rotateJacobi wrapper transformation is the same as using a Jacobi matrix as transformation matrix
    //  Specify some Jacobi rotation parameters
    size_t p = 2;
    size_t q = 5;
    double theta = 56.81;

    //  Use a Jacobi rotation matrix as transformation matrix
    Eigen::MatrixXd J = libwint::transformations::jacobiRotationMatrix(p, q, theta, K);

    Eigen::MatrixXd h = so_basis.get_h_SO();
    Eigen::Tensor<double, 4> g = so_basis.get_g_SO();

    Eigen::MatrixXd h_transformed_by_jacobi_matrix = libwint::transformations::transformOneElectronIntegrals(h, J);
    Eigen::Tensor<double, 4> g_transformed_by_jacobi_matrix = libwint::transformations::transformTwoElectronIntegrals(g, J);

    //  Rotate the SO basis
    so_basis.rotateJacobi(p, q, theta);


    BOOST_REQUIRE(h_transformed_by_jacobi_matrix.isApprox(so_basis.get_h_SO(), 1.0e-12));
    BOOST_REQUIRE(cpputil::linalg::areEqual(g_transformed_by_jacobi_matrix, so_basis.get_g_SO(), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( fcidump_constructor ) {

    libwint::SOBasis so_basis ("../tests/ref_data/beh_cation_631g_caitlin.FCIDUMP", 16);

    // Check if the one-electron integrals are read in correctly from a previous implementation
    Eigen::MatrixXd h_SO = so_basis.get_h_SO();

    BOOST_CHECK(std::abs(h_SO(0,0) - (-8.34082)) < 1.0e-5);
    BOOST_CHECK(std::abs(h_SO(5,1) - 0.381418) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(14,0) - 0.163205) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(13,6) - (-5.53204e-16)) < 1.0e-16);
    BOOST_CHECK(std::abs(h_SO(15,11) - (-0.110721)) < 1.0e-6);


    // Check if the two-electron integrals are read in correctly from a previous implementation
    Eigen::Tensor<double, 4> g_SO = so_basis.get_g_SO();

    BOOST_CHECK(std::abs(g_SO(2,5,4,4) - 0.0139645) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(2,6,3,0) - 5.16622e-18) < 1.0e-17);
    BOOST_CHECK(std::abs(g_SO(3,1,3,0) - (-0.0141251)) <  1.0e-6);
    BOOST_CHECK(std::abs(g_SO(4,6,4,6) - 0.0107791) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(4,15,11,1) - (9.33375e-19)) < 1.0e-17);
    BOOST_CHECK(std::abs(g_SO(6,10,5,9) - (-3.81422e-18)) < 1.0e-17);
    BOOST_CHECK(std::abs(g_SO(7,7,2,1) - (-0.031278)) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(8,15,9,9) - (-2.80093e-17)) < 1.0e-16);
    BOOST_CHECK(std::abs(g_SO(9,14,0,9) - 0.00161985) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(10,1,4,3) - 0.00264603) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(11,4,9,3) - (-0.0256623)) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(12,9,0,4) - 0.0055472) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(13,15,15,13) - 0.00766898) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(14,2,12,3) - 0.0104266) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(15,5,10,10) - 0.00562608) < 1.0e-7);
}


BOOST_AUTO_TEST_CASE ( fcidump_constructor_horton ) {

    // Check the same reference value as horton does
    libwint::SOBasis so_basis ("../tests/ref_data/h2_psi4_horton.FCIDUMP", 10);
    Eigen::Tensor<double, 4> g_SO = so_basis.get_g_SO();
    BOOST_CHECK(std::abs(g_SO(6,5,1,0) - 0.0533584656) <  1.0e-7);
}
