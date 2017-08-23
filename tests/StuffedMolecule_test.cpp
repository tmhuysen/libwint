#define BOOST_ALL_DYN_LINK

#define BOOST_TEST_MODULE "StuffedMolecule class"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain
#include "StuffedMolecule.hpp"


BOOST_AUTO_TEST_CASE( constructor ){

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    Molecule water (xyzfilename);
    BOOST_CHECK_EQUAL(water.natoms(), 3);


}