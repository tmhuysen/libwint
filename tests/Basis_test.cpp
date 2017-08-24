#define BOOST_ALL_DYN_LINK

#define BOOST_TEST_MODULE "Molecule class"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain
#include "Molecule.hpp"
#include "Basis.hpp"


BOOST_AUTO_TEST_CASE( constructor ){

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    std::string basis_name = "STO-3G";
    Molecule water (xyzfilename);

    Basis basis (water, basis_name);

    BOOST_CHECK_EQUAL(basis.molecule.natoms(), 3);
    BOOST_CHECK_EQUAL(basis.name, "STO-3G");
}