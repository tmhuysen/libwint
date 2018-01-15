#define BOOST_TEST_MODULE "transformations"

#include "Basis.hpp"
#include "transformations.hpp"
#include "utility.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE ( transform_one_electron_trivial ) {

    // Let's test a trivial transformation: i.e. with C being a unit matrix
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);

    Eigen::MatrixXd h_transformed = libwint::transform_one_electron_integrals(h, T);

    BOOST_CHECK(h_transformed.isApprox(h, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transform_one_electron_and_back ) {

    // Let's test if, if we transform h to h_transformed and back to h, we get the same result
    // We'll do this with random matrices
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    Eigen::MatrixXd T = Eigen::MatrixXd::Random(3, 3);  // The probability of a random matrix being singular is approximately 0
    Eigen::MatrixXd T_inverse = T.inverse();

    // Transform to SO basis
    Eigen::MatrixXd h_transformed = libwint::transform_one_electron_integrals(h, T);

    BOOST_CHECK(h.isApprox(libwint::transform_one_electron_integrals(h_transformed, T_inverse), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transform_two_electron_trivial ) {

    // Let's test a trivial transformation: i.e. with T being a unit matrix
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    Eigen::Tensor<double, 4> g (3, 3, 3, 3);
    g.setRandom();

    Eigen::Tensor<double, 4> g_transformed = libwint::transform_AO_to_SO(g, T);

    BOOST_CHECK(libwint::utility::are_equal(g, g_transformed, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transform_two_electron_olsens ) {

    // We can find a reference algorithm in the olsens branch from Ayer's lab
    Eigen::Tensor<double, 4> g_transformed_ref (2, 2, 2, 2);
    libwint::utility::read_array_from_file("../tests/ref_data/rotated1.data", g_transformed_ref);

    // Set the same matrix and tensor
    Eigen::MatrixXd T (2, 2);
    T << 1, 2, 3, 4;

    Eigen::Tensor<double, 4> g_transformed (2, 2, 2, 2);
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    g_transformed(i, j, k, l) = l + 2*k + 4*j + 8*i;
                }
            }
        }
    }

    Eigen::Tensor<double, 4> g_SO = libwint::transform_AO_to_SO(g_transformed, T);
    BOOST_CHECK(libwint::utility::are_equal(g_SO, g_transformed_ref, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( rotate ) {

    // Since the previous tests check the workings of the base function, we can check the workings of the wrapper function separately

    // Test if the functions don't accept non-unitary matrices
    Eigen::MatrixXd T (2, 2);
    T << 1, 2, 3, 4;

    Eigen::MatrixXd h = Eigen::MatrixXd::Random(2, 2);
    Eigen::Tensor<double, 4> g (2, 2, 2, 2);
    g.setRandom();

    BOOST_REQUIRE_THROW(libwint::rotate_integrals(h, T), std::invalid_argument);
    BOOST_REQUIRE_THROW(libwint::rotate_integrals(g, T), std::invalid_argument);

    // Test if the functions accept unitary matrices
    Eigen::MatrixXd U = Eigen::MatrixXd::Identity (2, 2);

    BOOST_REQUIRE_NO_THROW(libwint::rotate_integrals(h, U));
    BOOST_REQUIRE_NO_THROW(libwint::rotate_integrals(g, U));
}


BOOST_AUTO_TEST_CASE ( jacobi_rotation_matrix ) {

    // We can't create a Jacobi matrix for P > Q
    BOOST_REQUIRE_THROW(libwint::jacobi_rotation_matrix(3, 2, 1.0, 5), std::invalid_argument);

    // P+1 and Q+1 should both be smaller than M
    BOOST_REQUIRE_THROW(libwint::jacobi_rotation_matrix(1, 5, 1.0, 4), std::invalid_argument);
    BOOST_REQUIRE_NO_THROW(libwint::jacobi_rotation_matrix(2, 3, 1.0, 4));

    // A random Jacobi matrix is unitary
    BOOST_CHECK(libwint::jacobi_rotation_matrix(4, 7, 6.9921, 10).isUnitary());
    BOOST_CHECK(libwint::jacobi_rotation_matrix(1, 9, 78.00166, 22).isUnitary());

    // Let's test the easiest Jacobi matrix, one with theta = pi/2 and dimension 2
    Eigen::MatrixXd J = libwint::jacobi_rotation_matrix(0, 1, boost::math::constants::half_pi<double>(), 2);
    BOOST_CHECK(std::abs(J(0,0) - 0) < 1.0e-12);
    BOOST_CHECK(std::abs(J(0,1) - 1) < 1.0e-12);
    BOOST_CHECK(std::abs(J(1,0) - (-1)) < 1.0e-12);
    BOOST_CHECK(std::abs(J(0,0) - 0) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_jacobi_transformations ) {

    // Using a Jacobi matrix as a transformation matrix, check the transformed integrals with the results from olsens

    // Get the initial one- and two-electron integrals from olsens
    Eigen::MatrixXd h_SO (6, 6);
    Eigen::Tensor<double, 4> g_SO (6, 6, 6, 6);

    libwrp::utility::read_array_from_file("../tests/ref_data/lih_hf_sto6g_oneint.data", h_SO);
    libwrp::utility::read_array_from_file("../tests/ref_data/lih_hf_sto6g_twoint.data", g_SO);



    // Specify some Jacobi parameters to test a Jacobi rotation
    size_t P = 2;
    size_t Q = 4;
    double theta = 56.71;
    Eigen::MatrixXd U = libwrp::jacobi_rotation_matrix(P, Q, theta, 6);



    // Are the rotated one-electron integrals the same?
    Eigen::MatrixXd h_SO_rotated = libwrp::transform_one_electron_integrals(h_SO, U);

    Eigen::MatrixXd h_SO_rotated_olsens (6, 6);
    libwrp::utility::read_array_from_file("../tests/ref_data/lih_hf_sto6g_oneint_rotated.data", h_SO_rotated_olsens);

    BOOST_CHECK(h_SO_rotated.isApprox(h_SO_rotated_olsens, 1.0e-6));


    // Are the rotated two-electron integrals the same?
    Eigen::Tensor<double, 4> g_SO_rotated = libwrp::transform_two_electron_integrals(g_SO, U);

    Eigen::Tensor<double, 4> g_SO_rotated_olsens (6, 6, 6, 6);
    libwrp::utility::read_array_from_file("../tests/ref_data/lih_hf_sto6g_twoints_rotated.data", g_SO_rotated_olsens);

    BOOST_CHECK(libwrp::utility::are_equal(g_SO_rotated, g_SO_rotated_olsens, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( analytical_jacobi_one_electron_toy ) {

    // Set a toy one-electron matrix (which should be self-adjoint)
    Eigen::MatrixXd h (2, 2);
    h << 1, 2, 2, 4;

    // Set the test Jacobi parameters and test with the analytical formula
    size_t P = 0;
    size_t Q = 1;
    double theta = boost::math::constants::half_pi<double>();
    Eigen::MatrixXd U = libwint::jacobi_rotation_matrix(P, Q, theta, 2);
    BOOST_CHECK(libwint::rotate_integrals(h, U).isApprox(libwint::rotate_one_electron_integrals_jacobi(h, P, Q, theta)));

    // Test with another angle theta
    theta = 7.09214;
    U = libwint::jacobi_rotation_matrix(P, Q, theta, 2);
    BOOST_CHECK(libwint::rotate_integrals(h, U).isApprox(libwint::rotate_one_electron_integrals_jacobi(h, P, Q, theta)));
}


BOOST_AUTO_TEST_CASE ( analytical_jacobi_one_electron_h2o ) {

    // Calculate the kinetic energy integrals for H2O@STO-3G
    libint2::initialize();
    const std::string xyzfilename = "../tests/ref_data/h2o.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    const std::string basis_name = "STO-3G";
    libwint::Molecule water (xyzfilename);
    libwint::Basis basis (water, basis_name);
    basis.compute_kinetic_integrals();
    Eigen::MatrixXd T = basis.T;
    libint2::finalize();

    // Set the test Jacobi parameters and test with the analytical formula
    size_t P = 3;
    size_t Q = 6;
    double theta = 62.7219;
    Eigen::MatrixXd U = libwint::jacobi_rotation_matrix(P, Q, theta, basis.nbf());
    BOOST_CHECK(libwint::rotate_integrals(T, U).isApprox(libwint::rotate_one_electron_integrals_jacobi(T, P, Q, theta)));
}
