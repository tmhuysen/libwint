#define BOOST_TEST_MODULE "Molecule"


#include "Molecule.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain



BOOST_AUTO_TEST_CASE ( parse_xyz_filename ) {

    // Make sure we get an error when a nonsense path is given for the .xyz file name
    BOOST_REQUIRE_THROW(libwint::Molecule ("this is a nonsense data path"), std::runtime_error);

    // Make sure we get an error when a path with a wrong extension is given
    BOOST_REQUIRE_THROW(libwint::Molecule ("../tests/ref_data/nuclear.data"), std::runtime_error);

    // Make sure we don't get an error when a correct path is given
    BOOST_REQUIRE_NO_THROW(libwint::Molecule ("../tests/ref_data/h2o.xyz"));
}


BOOST_AUTO_TEST_CASE ( distance ) {

    // Create a fictitious molecule from some libint2::Atoms (charge, x, y ,z)
    libint2::Atom A {0, 0, 3, 0};
    libint2::Atom B {0, 0, 0, 4};
    libint2::Atom C {0, 3, 0, 0};
    libint2::Atom D {0, 0, 0, 5};
    std::vector<libint2::Atom> atoms = {A, B, C, D};
    libwint::Molecule molecule (atoms);


    // Check their distances
    BOOST_CHECK(std::abs(molecule.calculateInternuclearDistance(0, 1) - 5) < 1.0e-12);
    BOOST_CHECK(std::abs(molecule.calculateInternuclearDistance(0, 2) - std::sqrt(18.0)) < 1.0e-12);
    BOOST_CHECK(std::abs(molecule.calculateInternuclearDistance(0, 1) - molecule.calculateInternuclearDistance(1, 2)) < 1.0e-12);
    BOOST_CHECK(std::abs(molecule.calculateInternuclearDistance(1, 2) - 5) < 1.0e-12);
    BOOST_CHECK(std::abs(molecule.calculateInternuclearDistance(1, 3) - 1) < 1.0e-12);

    // Check that the distances are symmetric
    BOOST_CHECK(std::abs(molecule.calculateInternuclearDistance(0, 1) - molecule.calculateInternuclearDistance(1, 0)) < 1.0e-12);
    BOOST_CHECK(std::abs(molecule.calculateInternuclearDistance(0, 2) - molecule.calculateInternuclearDistance(2, 0)) < 1.0e-12);
    BOOST_CHECK(std::abs(molecule.calculateInternuclearDistance(1, 2) - molecule.calculateInternuclearDistance(2, 1)) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( molecule_constructor ) {

    // Create some molecule objects
    const std::string xyzfilename = "../tests/ref_data/h2o.xyz";  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::Molecule water (xyzfilename);
    libwint::Molecule water_anion (xyzfilename, -1);
    libwint::Molecule water_neutral (xyzfilename, 0);
    libwint::Molecule water_cation (xyzfilename, +1);

    // Test the number of electrons created by the constructor
    BOOST_CHECK_EQUAL(water.get_N(), 10);
    BOOST_CHECK_EQUAL(water_anion.get_N(), 11);
    BOOST_CHECK_EQUAL(water_neutral.get_N(), 10);
    BOOST_CHECK_EQUAL(water_cation.get_N(), 9);
}


BOOST_AUTO_TEST_CASE ( methods_water ) {

    // We have reference internuclear repulsion energy from Horton
    double ref_internuclear_repulsion_energy = 8.00236693455;

    // Create the water molecule
    libwint::Molecule water ("../tests/ref_data/h2o.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory

    // Test the basic methods
    BOOST_CHECK_EQUAL(water.numberOfAtoms(), 3);
    BOOST_CHECK_EQUAL(water.calculateTotalNucleicCharge(), 10);

    // Test the calculation of the nuclear repulsion energy
    BOOST_CHECK(std::abs(water.calculateInternuclearRepulsionEnergy() - ref_internuclear_repulsion_energy) < 1.0e-07);  // reference data from horton
}


BOOST_AUTO_TEST_CASE ( methods_h2 ) {

    // We have reference internuclear repulsion energy from Horton
    double ref_internuclear_repulsion_energy = 0.714285658963;

    // Create the hydrogen gas molecule
    libwint::Molecule h2 ("../tests/ref_data/h2.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory

    // Test the basic methods
    BOOST_CHECK_EQUAL(h2.numberOfAtoms(), 2);
    BOOST_CHECK_EQUAL(h2.calculateTotalNucleicCharge(), 2);

    // Test the calculation of the nuclear repulsion energy
    BOOST_CHECK(std::abs(h2.calculateInternuclearRepulsionEnergy() - ref_internuclear_repulsion_energy) < 1.0e-07);  // reference data from horton
}
