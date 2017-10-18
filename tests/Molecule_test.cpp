#define BOOST_TEST_MODULE "Molecule"

#include "Molecule.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE( constructor ) {

    // Create the molecules
    constexpr auto xyzfilename = "../../docs/h2o.xyz"; // Anticipate an out-of source build, so we need one level higher in directories
    Wrapper::Molecule water (xyzfilename);
    Wrapper::Molecule water_anion (xyzfilename, -1);
    Wrapper::Molecule water_neutral (xyzfilename, 0);
    Wrapper::Molecule water_cation (xyzfilename, +1);

    // Test the constructor
    BOOST_CHECK_EQUAL(water.xyz_filename, "../../docs/h2o.xyz");

    // Test the number of electrons created by the constructor
    BOOST_CHECK_EQUAL(water.nelec, 10);
    BOOST_CHECK_EQUAL(water_anion.nelec, 11);
    BOOST_CHECK_EQUAL(water_neutral.nelec, 10);
    BOOST_CHECK_EQUAL(water_cation.nelec, 9);
}


BOOST_AUTO_TEST_CASE ( methods_water ) {

    // Create the water molecule
    constexpr auto xyzfilename = "../../docs/h2o.xyz"; // Anticipate an out-of source build, so we need one level higher in directories
    Wrapper::Molecule water (xyzfilename);

    // Test the basic methods
    BOOST_CHECK_EQUAL(water.natoms(), 3);
    BOOST_CHECK_EQUAL(water.nucleic_charge(), 10);

    // Test the calculation of the nuclear repulsion energy
    BOOST_CHECK(std::abs(water.internuclear_repulsion() - 8.00236693455) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( methods_h2 ) {

    // Create the water molecule
    constexpr auto xyzfilename = "../../docs/h2.xyz"; // Anticipate an out-of source build, so we need one level higher in directories
    Wrapper::Molecule h2 (xyzfilename);

    // Test the basic methods
    BOOST_CHECK_EQUAL(h2.natoms(), 3);
    BOOST_CHECK_EQUAL(h2.nucleic_charge(), 2);

    // Test the calculation of the nuclear repulsion energy
    BOOST_CHECK(std::abs(h2.internuclear_repulsion() - 0.714285658963) < 1.0e-08);
}