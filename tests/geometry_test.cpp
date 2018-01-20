#define BOOST_TEST_MODULE "geometry"

#include "geometry.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE ( distance ) {

    // Create some libint2::Atoms (charge, x, y ,z)
    libint2::Atom A {0, 0, 3, 0};
    libint2::Atom B {0, 0, 0, 4};
    libint2::Atom C {0, 3, 0, 0};
    libint2::Atom D {0, 0, 0, 5};

    // Check their distances
    BOOST_CHECK(std::abs(libwint::distance(A, B) - 5) < 1.0e-8);
    BOOST_CHECK(std::abs(libwint::distance(A, C) - std::sqrt(18.0)) < 1.0e-8);
    BOOST_CHECK(std::abs(libwint::distance(A, B) - libwint::distance(B, C)) < 1.0e-8);
    BOOST_CHECK(std::abs(libwint::distance(B, C) - 5) < 1.0e-8);
    BOOST_CHECK(std::abs(libwint::distance(B, D) - 1) < 1.0e-8);

    // Check that the distances are symmetric
    BOOST_CHECK(std::abs(libwint::distance(A, B) - libwint::distance(B, A)) < 1.0e-8);
    BOOST_CHECK(std::abs(libwint::distance(A, C) - libwint::distance(C, A)) < 1.0e-8);
    BOOST_CHECK(std::abs(libwint::distance(B, C) - libwint::distance(C, B)) < 1.0e-8);
}