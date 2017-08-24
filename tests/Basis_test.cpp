#define BOOST_ALL_DYN_LINK

#define BOOST_TEST_MODULE "Molecule class"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain
#include "Molecule.hpp"
#include "Basis.hpp"
#include <libint2.hpp>


BOOST_AUTO_TEST_CASE( constructor ){
    // Initialize libint2
    libint2::initialize();

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    std::string basis_name = "STO-3G";
    Molecule water (xyzfilename);

    Basis basis (water, basis_name);


    BOOST_CHECK_EQUAL(basis.molecule.natoms(), 3);
    BOOST_CHECK_EQUAL(basis.name, "STO-3G");

    // Finalize libint2
    libint2::finalize();
}

BOOST_AUTO_TEST_CASE( integrals ){
    // Initialize libint2
    libint2::initialize();

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    std::string basis_name = "STO-3G";
    Molecule water (xyzfilename);

    Basis basis (water, basis_name);

    auto S = basis.compute_overlap_integrals();
    auto T = basis.compute_kinetic_integrals();
    auto V = basis.compute_nuclear_integrals();

    auto tei = basis.compute_two_electron_integrals();

    // Finalize libint2
    libint2::finalize();
}


