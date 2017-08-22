//
// Created by Laurent Lemmens on 21/08/17.
//

#define BOOST_TEST_MODULE "SOME TEST GOING ON"
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain
#include <libint2.hpp>

BOOST_AUTO_TEST_CASE( my_test )
{

    // 1. MOLECULE & BASIS SET SPECIFICATION
    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    std::ifstream input_file(xyzfilename);
    auto atoms = libint2::read_dotxyz(input_file);
    libint2::BasisSet obs ("STO-3G", atoms);    // obs: orbital basis set
                                                // a libint2::BasisSet is a decorated std::vector<libint2::Shell>

    BOOST_CHECK(obs.size() == 5);
    BOOST_CHECK(obs.nbf() == 7);

}