#define BOOST_TEST_MODULE "transformations"

#include "transformations.hpp"
#include "utility.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE ( one_electron ) {

    // Since I don't have any data (yet) that can confirm the AO to SO integral transformations, we should at least check that the code works with random matrices

    Eigen::MatrixXd C = Eigen::MatrixXd::Random(3, 3);
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);

    BOOST_REQUIRE_NO_THROW(libwrp::transform_AO_integrals_to_SO(h, C));
}


BOOST_AUTO_TEST_CASE ( two_electron_trivial ) {

    // Let's test a trivial transformation: i.e. with C being a unit matrix
    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(3, 3);
    Eigen::Tensor<double, 4> g_AO (3, 3, 3, 3);
    g_AO.setRandom();

    Eigen::Tensor<double, 4> g_SO = libwrp::transform_AO_integrals_to_SO(g_AO, C);

    BOOST_CHECK(libwrp::utility::are_equal(g_AO, g_SO, 1.0e-12));
}

BOOST_AUTO_TEST_CASE ( two_electron_olsens1 ) {

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

    Eigen::Tensor<double, 4> g_SO = libwrp::transform_AO_integrals_to_SO(g_AO, C);

    libwrp::utility::print(rotated_tei_ref);
    libwrp::utility::print(g_SO);

    BOOST_CHECK(libwrp::utility::are_equal(g_SO, rotated_tei_ref, 1.0e-06));
}
