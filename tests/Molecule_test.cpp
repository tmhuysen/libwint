#define BOOST_TEST_MODULE "Molecule"

#include "Molecule.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE( constructor ){

    constexpr auto xyzfilename = "../../docs/h2o.xyz"; // Anticipate an out-of source build, so we need one level higher in directories
    Wrapper::Molecule water (xyzfilename);

    BOOST_CHECK_EQUAL(water.xyz_filename, "../../docs/h2o.xyz");
    BOOST_CHECK_EQUAL(water.natoms(), 3);

    BOOST_CHECK_EQUAL(water.nucleic_charge(), 10);
    BOOST_CHECK_EQUAL(water.nelec, 10);
}