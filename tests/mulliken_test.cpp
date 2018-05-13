#define BOOST_TEST_MODULE "SOMullikenBasis_test"

#include "SOMullikenBasis.hpp"
#include "transformations.hpp"

#include <cpputil.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( mulliken_test ) {
    Eigen::Matrix<double,7,7> aa;
    aa << 0.999999, 2.75445e-05, 2.10116e-16, -8.05074e-06, -2.79652e-18, 4.63697e-05, 3.64526e-17,
        2.75445e-05, 0.995582, -2.86904e-17, -0.00753316, 8.17833e-17,  0.00425522, -1.91957e-16,
        2.10116e-16, -2.86904e-17, 0.976792, 2.28438e-17 , 3.71927e-17,  3.65174e-16, -0.00047531,
        -8.05074e-06, -0.00753316, 2.28438e-17, 0.980476, 6.4507e-17, -0.0182836, 1.26353e-15,
        -2.79652e-18, 8.17833e-17, 3.71927e-17, 6.4507e-17, 0.999083, -4.46521e-15, -2.88329e-16,
        4.63697e-05, 0.00425522,  3.65174e-16, -0.0182836, -4.46521e-15,  0.0247316, -3.84334e-17,
        3.64526e-17, -1.91957e-16, -0.00047531, 1.26353e-15, -2.88329e-16, -3.84334e-17 , 0.0233358;
    Eigen::Matrix<double,7,7> C;
    C <<  -0.994435, -0.239159, 2.70425e-16, 0.0936834, 3.0127e-32, -0.11164, 1.19708e-15,
            -0.024097, 0.885736, -1.51665e-15, -0.479587, -4.1368e-31, 0.669578, -7.13724e-15,
            -7.16711e-19, -1.69007e-16, 0.607285, -2.60703e-15, 1.15279e-15, 9.63173e-15, 0.919233,
            -0.00316155, 0.0858963, 1.81567e-15, 0.74743, 2.31276e-31, 0.73849, -6.72651e-15,
            -2.04115e-34, 8.53707e-32, 1.68919e-15, -5.09334e-30, -1, -1.59442e-30, -1.14861e-16,
            0.00459374, 0.144039, 0.452998, 0.329472, 6.92404e-16, -0.709849, -0.732461,
            0.00459374, 0.144039, -0.452998, 0.329472, -6.92404e-16, -0.709849, 0.732461;

    Eigen::MatrixXd bb = aa;
    libwint::Molecule water ("../tests/ref_data/h2o.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();
    libwint::SOMullikenBasis so_basis (ao_basis, C);
    so_basis.calculateMullikenMatrix({0,1});

    // TO:DO add ref
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE ( transform_jacobi_mulliken ) {

    // Create an SOBasis instance with a coefficient matrix being the identity matrix (little hack that we can use to test transformations)
    libwint::Molecule water ("../tests/ref_data/h2o.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();
    size_t K = ao_basis.calculateNumberOfBasisFunctions();

    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(K, K);
    libwint::SOMullikenBasis so_basis (ao_basis, C);
    libwint::SOMullikenBasis so_basis_check (ao_basis, C);
    so_basis.calculateMullikenMatrix({0,1});

    // Test if the rotateJacobi wrapper transformation is the same as using a Jacobi matrix as transformation matrix
    //  Specify some Jacobi rotation parameters
    size_t p = 2;
    size_t q = 5;
    double theta = 56.81;

    //  Use a Jacobi rotation matrix as transformation matrix
    Eigen::MatrixXd J = libwint::transformations::jacobiRotationMatrix(p, q, theta, K);

    Eigen::MatrixXd h = so_basis.get_h_SO();
    Eigen::MatrixXd m = so_basis.get_mulliken_matrix();
    Eigen::Tensor<double, 4> g = so_basis.get_g_SO();

    Eigen::MatrixXd h_transformed_by_jacobi_matrix = libwint::transformations::transformOneElectronIntegrals(h, J);
    Eigen::MatrixXd m_transformed_by_jacobi_matrix = libwint::transformations::transformOneElectronIntegrals(m, J);
    Eigen::Tensor<double, 4> g_transformed_by_jacobi_matrix = libwint::transformations::transformTwoElectronIntegrals(g, J);

    //  Rotate the SO basis
    so_basis.rotateJacobi(p, q, theta);


    BOOST_REQUIRE(h_transformed_by_jacobi_matrix.isApprox(so_basis.get_h_SO(), 1.0e-6));
    BOOST_REQUIRE(m_transformed_by_jacobi_matrix.isApprox(so_basis.get_mulliken_matrix(), 1.0e-6));
    BOOST_REQUIRE(cpputil::linalg::areEqual(g_transformed_by_jacobi_matrix, so_basis.get_g_SO(), 1.0e-6));
}

BOOST_AUTO_TEST_CASE ( mulliken_copy ) {

    // Create an SOBasis instance with a coefficient matrix being the identity matrix (little hack that we can use to test transformations)
    libwint::Molecule water ("../tests/ref_data/h2o.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();
    size_t K = ao_basis.calculateNumberOfBasisFunctions();

    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(K, K);
    libwint::SOMullikenBasis so_basis (ao_basis, C);
    libwint::SOMullikenBasis so_basis_check (ao_basis, C);
    so_basis.calculateMullikenMatrix({0,1});
    so_basis.rotateJacobi(1,2,30);
    // Make sure we don't get an error when copying
    BOOST_REQUIRE_NO_THROW(so_basis_check.copy(so_basis));

    BOOST_CHECK(so_basis_check.get_mulliken_matrix().isApprox(so_basis.get_mulliken_matrix(), 1.0e-6));


}



BOOST_AUTO_TEST_CASE ( reader ) {

    // Create an SOBasis instance with a coefficient matrix being the identity matrix (little hack that we can use to test transformations)
    //  Start reading in the one- and two-electron integrals
    std::string fcidump_filename = "../tests/ref_data/no_0.5_PB";
    libwint::SOMullikenBasis so_basis(fcidump_filename,10);
    BOOST_CHECK(true);


}