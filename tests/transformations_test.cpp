#define BOOST_TEST_MODULE "transformations"

#include "transformations.hpp"
#include "utility.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE ( transform_one_electron_trivial ) {

    // Let's test a trivial transformation: i.e. with C being a unit matrix
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);

    Eigen::MatrixXd h_transformed = libwrp::transform_one_electron_integrals(h, T);

    BOOST_CHECK(h_transformed.isApprox(h, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transform_one_electron_and_back ) {

    // Let's test if, if we transform h to h_transformed and back to h, we get the same result
    // We'll do this with random matrices
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    Eigen::MatrixXd T = Eigen::MatrixXd::Random(3, 3);  // The probability of a random matrix being singular is approximately 0
    Eigen::MatrixXd T_inverse = T.inverse();

    // Transform to SO basis
    Eigen::MatrixXd h_transformed = libwrp::transform_one_electron_integrals(h, T);

    BOOST_CHECK(h.isApprox(libwrp::transform_one_electron_integrals(h_transformed, T_inverse), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transform_two_electron_trivial ) {

    // Let's test a trivial transformation: i.e. with T being a unit matrix
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    Eigen::Tensor<double, 4> g (3, 3, 3, 3);
    g.setRandom();

    Eigen::Tensor<double, 4> g_transformed = libwrp::transform_AO_to_SO(g, T);

    BOOST_CHECK(libwrp::utility::are_equal(g, g_transformed, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transform_two_electron_olsens ) {

    // We can find a reference algorithm in the olsens branch from Ayer's lab
    Eigen::Tensor<double, 4> g_transformed_ref (2, 2, 2, 2);
    libwrp::utility::read_array_from_file("../tests/ref_data/rotated1.data", g_transformed_ref);

    // Set the same matrix and tensor
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

    Eigen::Tensor<double, 4> g_SO = libwrp::transform_AO_to_SO(g_transformed, T);
    BOOST_CHECK(libwrp::utility::are_equal(g_SO, g_transformed_ref, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( rotate ) {

    // Since the previous tests check the workings of the base function, we can check the workings of the wrapper function separately

    // Test if the functions don't accept non-unitary matrices
    Eigen::MatrixXd T (2, 2);
    T << 1, 2, 3, 4;

    Eigen::MatrixXd h = Eigen::MatrixXd::Random(2, 2);
    Eigen::Tensor<double, 4> g (2, 2, 2, 2);
    g.setRandom();

    BOOST_REQUIRE_THROW(libwrp::rotate_integrals(h, T), std::invalid_argument);
    BOOST_REQUIRE_THROW(libwrp::rotate_integrals(g, T), std::invalid_argument);

    // Test if the functions accept unitary matrices
    Eigen::MatrixXd U = Eigen::MatrixXd::Identity (2, 2);

    BOOST_REQUIRE_NO_THROW(libwrp::rotate_integrals(h, U));
    BOOST_REQUIRE_NO_THROW(libwrp::rotate_integrals(g, U));
}


BOOST_AUTO_TEST_CASE ( jacobi_rotation_matrix ) {

    // We can't create a Jacobi matrix for P > Q
    BOOST_REQUIRE_THROW(libwrp::jacobi_rotation_matrix(3, 2, 1.0, 5), std::invalid_argument);

    // P+1 and Q+1 should both be smaller than M
    BOOST_REQUIRE_THROW(libwrp::jacobi_rotation_matrix(1, 5, 1.0, 4), std::invalid_argument);
    BOOST_REQUIRE_NO_THROW(libwrp::jacobi_rotation_matrix(2, 3, 1.0, 4));

    // A random Jacobi matrix is unitary
    BOOST_CHECK(libwrp::jacobi_rotation_matrix(4, 7, 6.9921, 10).isUnitary());
    BOOST_CHECK(libwrp::jacobi_rotation_matrix(1, 9, 78.00166, 22).isUnitary());

}
