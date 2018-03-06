#define BOOST_TEST_MODULE "SOBasis"


#include "SOBasis.hpp"

#include "transformations.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain



BOOST_AUTO_TEST_CASE ( transform_jacobi ) {

    // Create an SOBasis instance with a coefficient matrix being the identity matrix (little hack that we can use to test transformations)
    libwint::Molecule water ("../tests/ref_data/h2o.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::AOBasis ao_basis (water, "STO-3G");
    size_t K = ao_basis.calculateNumberOfBasisFunctions();
    ao_basis.calculateIntegrals();

    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(K, K);
    libwint::SOBasis so_basis (ao_basis, C);


    // Test if the transformJacobi wrapper transformation is the same as using a Jacobi matrix as transformation matrix
    //  Specify some Jacobi rotation parameters
    size_t p = 2;
    size_t q = 5;
    double theta = 56.81;

    //  Use a Jacobi rotation matrix as transformation matrix
    Eigen::MatrixXd J = libwint::transformations::jacobiRotationMatrix(p, q, theta, K);

    Eigen::MatrixXd h = so_basis.get_h_SO();
    Eigen::Tensor<double, 4> g = so_basis.get_g_SO();

    Eigen::MatrixXd h_transformed_by_jacobi_matrix = libwint::transformations::transformOneElectronIntegrals(h, J);
    Eigen::Tensor<double, 4> g_transformed_by_jacobi_matrix = libwint::transformations::transformTwoElectronIntegrals(g, J);

    //  Rotate the SO basis
    so_basis.transformJacobi(p, q, theta);


    BOOST_REQUIRE(h_transformed_by_jacobi_matrix.isApprox(so_basis.get_h_SO(), 1.0e-12));
    BOOST_REQUIRE(cpputil::linalg::areEqual(g_transformed_by_jacobi_matrix, so_basis.get_g_SO(), 1.0e-12));
}
