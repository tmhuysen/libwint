#define BOOST_TEST_MODULE "transformations"

#include "transformations.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE ( one_electron ) {

    // Since I don't have any data (yet) that can confirm the AO to SO integral transformations, we should at least check that the code works with random matrices
    
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(3, 3);
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);

    BOOST_REQUIRE_NO_THROW(libwrp::transform_AO_integrals_to_SO(h, C));
}


BOOST_AUTO_TEST_CASE ( two_electron ) {

    // Since I don't have any data (yet) that can confirm the AO to SO integral transformations, we should at least check that the code works with random tensors and matrices

    Eigen::MatrixXd C = Eigen::MatrixXd::Random(3, 3);
    Eigen::Tensor<double, 4> tei (3, 3, 3, 3);
    tei.setRandom();

    BOOST_REQUIRE_NO_THROW(libwrp::transform_AO_integrals_to_SO(tei, C));
}