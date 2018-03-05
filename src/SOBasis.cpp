#include "SOBasis.hpp"

#include "transformations.hpp"



namespace libwint {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param atomic orbital instance and a coefficient matrix @param C (i.e. a basis transformation matrix) that links the SO basis to the AO basis
 */
SOBasis::SOBasis(const libwint::AOBasis& ao_basis, const Eigen::MatrixXd& C):
        K(ao_basis.calculateNumberOfBasisFunctions())
{
    Eigen::MatrixXd h_AO = ao_basis.get_T() + ao_basis.get_V();
    this->h_SO = libwint::transformations::transform_AO_to_SO(h_AO, C);

    this->g_SO = libwint::transformations::transform_AO_to_SO(ao_basis.get_g(), C);
}


/*
 *  PUBLIC METHODS
 */

/**
 *  Transform the one- and two-electron integrals according to the basis transformation matrix @param T
 */
void SOBasis::transform(const Eigen::MatrixXd& T) {

    this->h_SO = libwint::transformations::transformOneElectronIntegrals(this->h_SO, T);
    this->g_SO = libwint::transformations::transformTwoElectronIntegrals(this->g_SO, T);
}


/**
 *  Transform the one- and two-electron integrals according to the Jacobi rotation parameters p, q and a given angle.
 */
void SOBasis::transformJacobi(size_t p, size_t q, double angle) {

    // Construct a Jacobi rotation matrix and to the transformation with that matrix
    auto dim = static_cast<size_t>(this->h_SO.cols());  // cols() returns a long
    Eigen::MatrixXd jacobi_rotation_matrix = libwint::transformations::jacobiRotationMatrix(p, q, angle, dim);

    this->transform(jacobi_rotation_matrix);
}


}  // namespace libwint
