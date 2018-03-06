#define BOOST_TEST_MODULE "transformations"

#include "transformations.hpp"

#include <boost/math/constants/constants.hpp>
#include "cpputil.hpp"

#include "AOBasis.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain



BOOST_AUTO_TEST_CASE ( transform_one_electron_trivial ) {

    // Let's test a trivial transformation: i.e. with C being a unit matrix
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);

    Eigen::MatrixXd h_transformed = libwint::transformations::transformOneElectronIntegrals(h, T);

    BOOST_CHECK(h_transformed.isApprox(h, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transform_one_electron_and_back ) {

    // Let's test if, if we transform h to h_transformed and back to h, we get the same result
    // We'll do this with random matrices
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    Eigen::MatrixXd T = Eigen::MatrixXd::Random(3, 3);  // the probability of a random matrix being singular is approximately 0
    Eigen::MatrixXd T_inverse = T.inverse();

    // Transform to SO basis
    Eigen::MatrixXd h_transformed = libwint::transformations::transformOneElectronIntegrals(h, T);

    BOOST_CHECK(h.isApprox(libwint::transformations::transformOneElectronIntegrals(h_transformed, T_inverse), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transform_two_electron_trivial ) {

    // Let's test a trivial transformation: i.e. with T being a unit matrix
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    Eigen::Tensor<double, 4> g (3, 3, 3, 3);
    g.setRandom();

    Eigen::Tensor<double, 4> g_transformed = libwint::transformations::transform_AO_to_SO(g, T);

    BOOST_CHECK(cpputil::linalg::areEqual(g, g_transformed, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transform_two_electron_olsens ) {

    // We can find a reference algorithm in the olsens module from Ayer's lab
    Eigen::Tensor<double, 4> g_transformed_ref (2, 2, 2, 2);
    cpputil::io::readArrayFromFile("../tests/ref_data/rotated_olsens.data", g_transformed_ref);

    // Set an example transformation matrix and two-electron integrals tensor
    Eigen::MatrixXd T (2, 2);
    T << 1, 2, 3, 4;

    Eigen::Tensor<double, 4> g_transformed (2, 2, 2, 2);
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    g_transformed(i, j, k, l) = l + 2*k + 4*j + 8*i;
                }
            }
        }
    }

    Eigen::Tensor<double, 4> g_SO = libwint::transformations::transform_AO_to_SO(g_transformed, T);
    BOOST_CHECK(cpputil::linalg::areEqual(g_SO, g_transformed_ref, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( transform_so_to_ao ) {

    // Check that AO->SO->AO is an identity operation
    // We'll do this with random matrices
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    Eigen::MatrixXd T = Eigen::MatrixXd::Random(3, 3);  // the probability of a random matrix being singular is approximately 0

    // Transform to SO basis
    Eigen::MatrixXd h_SO = libwint::transformations::transform_AO_to_SO(h, T);

    BOOST_CHECK(h.isApprox(libwint::transformations::transform_SO_to_AO(h_SO, T), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( jacobi_rotation_matrix ) {

    // We can't create a Jacobi matrix for P > Q
    BOOST_REQUIRE_THROW(libwint::transformations::jacobiRotationMatrix(3, 2, 1.0, 5), std::invalid_argument);

    // P+1 and Q+1 should both be smaller than M
    BOOST_REQUIRE_THROW(libwint::transformations::jacobiRotationMatrix(1, 5, 1.0, 4), std::invalid_argument);
    BOOST_REQUIRE_NO_THROW(libwint::transformations::jacobiRotationMatrix(2, 3, 1.0, 4));

    // A random Jacobi matrix is unitary
    BOOST_CHECK(libwint::transformations::jacobiRotationMatrix(4, 7, 6.9921, 10).isUnitary());
    BOOST_CHECK(libwint::transformations::jacobiRotationMatrix(1, 9, 78.00166, 22).isUnitary());

    // Let's test the easiest Jacobi matrix, one with theta = pi/2 and dimension 2
    Eigen::MatrixXd J = libwint::transformations::jacobiRotationMatrix(0, 1, boost::math::constants::half_pi<double>(), 2);
    std::cout << J << std::endl;
    BOOST_CHECK(std::abs(J(0,0) - 0) < 1.0e-12);
    BOOST_CHECK(std::abs(J(0,1) - 1) < 1.0e-12);
    BOOST_CHECK(std::abs(J(1,0) - (-1)) < 1.0e-12);
    BOOST_CHECK(std::abs(J(0,0) - 0) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_jacobi_transformations ) {

    // Using a Jacobi matrix as a transformation matrix, check the transformed integrals with the results from olsens

    // Get the initial one- and two-electron integrals from olsens
    Eigen::MatrixXd h_SO (6, 6);
    Eigen::Tensor<double, 4> g_SO (6, 6, 6, 6);

    cpputil::io::readArrayFromFile("../tests/ref_data/lih_hf_sto6g_oneint.data", h_SO);
    cpputil::io::readArrayFromFile("../tests/ref_data/lih_hf_sto6g_twoint.data", g_SO);


    // Specify some Jacobi parameters to test a Jacobi rotation
    size_t p = 2;
    size_t q = 4;
    double theta = 56.71;

    // Are the rotated one-electron integrals the same?
    Eigen::MatrixXd h_SO_rotated = libwint::transformations::rotateOneElectronIntegralsJacobi(h_SO, p, q, theta);

    Eigen::MatrixXd h_SO_rotated_olsens (6, 6);
    cpputil::io::readArrayFromFile("../tests/ref_data/lih_hf_sto6g_oneint_rotated.data", h_SO_rotated_olsens);
    BOOST_CHECK(h_SO_rotated.isApprox(h_SO_rotated_olsens, 1.0e-6));


    // Are the rotated two-electron integrals the same?
    Eigen::Tensor<double, 4> g_SO_rotated = libwint::transformations::rotateTwoElectronIntegralsJacobi(g_SO, p, q, theta);

    Eigen::Tensor<double, 4> g_SO_rotated_olsens (6, 6, 6, 6);
    cpputil::io::readArrayFromFile("../tests/ref_data/lih_hf_sto6g_twoints_rotated.data", g_SO_rotated_olsens);
    BOOST_CHECK(cpputil::linalg::areEqual(g_SO_rotated, g_SO_rotated_olsens, 1.0e-06));
}
