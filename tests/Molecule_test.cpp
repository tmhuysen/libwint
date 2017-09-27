#define BOOST_ALL_DYN_LINK

#define BOOST_TEST_MODULE "Molecule"

#include "Molecule.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE( constructor ){

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    Wrapper::Molecule water (xyzfilename);

    BOOST_CHECK_EQUAL(water.xyz_filename, "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz");
    BOOST_CHECK_EQUAL(water.natoms(), 3);
}