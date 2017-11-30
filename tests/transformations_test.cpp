#define BOOST_TEST_MODULE "transformations"

#include "transformations.hpp"
#include "utility.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE ( one_electron_trivial ) {

    // Let's test a trivial transformation: i.e. with C being a unit matrix
    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd h_AO = Eigen::MatrixXd::Random(3, 3);

    Eigen::MatrixXd h_SO = libwrp::transform_AO_to_SO(h_AO, C);

    BOOST_CHECK(h_AO.isApprox(h_SO, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( one_electron_AO_SO_back ) {

    // Let's test if, if we transform f_AO to f_SO and back to f_AO, we get the same result
    // We'll do this with random matrices
    Eigen::MatrixXd f_AO = Eigen::MatrixXd::Random(3, 3);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(3, 3);  // The probability of a random matrix being singular is approximately 0

    // Transform to SO basis
    Eigen::MatrixXd f_SO = libwrp::transform_AO_to_SO(f_AO, C);

    BOOST_CHECK(f_AO.isApprox(libwrp::transform_SO_to_AO(f_SO, C), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( two_electron_trivial ) {

    // Let's test a trivial transformation: i.e. with C being a unit matrix
    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(3, 3);
    Eigen::Tensor<double, 4> g_AO (3, 3, 3, 3);
    g_AO.setRandom();

    Eigen::Tensor<double, 4> g_SO = libwrp::transform_AO_to_SO(g_AO, C);

    BOOST_CHECK(libwrp::utility::are_equal(g_AO, g_SO, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( two_electron_olsens ) {

    // We can find a reference algorithm in the olsens branch from Ayer's lab
    Eigen::Tensor<double, 4> rotated_tei_ref (2, 2, 2, 2);
    libwrp::utility::read_array_from_file("../tests/ref_data/rotated1.data", rotated_tei_ref);

    // Set the same matrix and tensor
    Eigen::MatrixXd C (2, 2);
    C << 1, 2, 3, 4;

    Eigen::Tensor<double, 4> g_AO (2, 2, 2, 2);
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    g_AO(i, j, k, l) = l + 2*k + 4*j + 8*i;
                }
            }
        }
    }

    Eigen::Tensor<double, 4> g_SO = libwrp::transform_AO_to_SO(g_AO, C);
    BOOST_CHECK(libwrp::utility::are_equal(g_SO, rotated_tei_ref, 1.0e-06));
}
