#define BOOST_ALL_DYN_LINK

#define BOOST_TEST_MODULE "Molecule class"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain
#include "Molecule.hpp"


BOOST_AUTO_TEST_CASE( constructor ){

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    Molecule water (xyzfilename);

    BOOST_CHECK_EQUAL(water.xyz_filename, "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz");
    BOOST_CHECK_EQUAL(water.natoms(), 3);
}