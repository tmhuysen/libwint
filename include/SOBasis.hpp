#ifndef LIBWINT_SOBASIS_HPP
#define LIBWINT_SOBASIS_HPP


#include <Eigen/Dense>

#include "AOBasis.hpp"



namespace libwint {


class SOBasis {
private:
    const size_t K;  // the number of spatial orbitals

    Eigen::MatrixXd h_SO;  // the one-electron integrals (core Hamiltonian) in the spatial orbital basis
    Eigen::Tensor<double, 4> g_SO;  // the two-electron repulsion integrals in the spatial orbital basis



public:
    // Constructors
    /**
     *  Constructor based on a given @param atomic orbital instance and a coefficient matrix @param C (i.e. a basis transformation matrix) that links the SO basis to the AO basis
     */
    SOBasis(const libwint::AOBasis& ao_basis, const Eigen::MatrixXd& C);


    // Getters
    Eigen::MatrixXd get_h_SO() const { return this->h_SO; }
    Eigen::Tensor<double, 4> get_g_SO() const { return this->g_SO; }
    const size_t get_K() const { return this->K; }


    // Methods
    /**
     *  Transform the one- and two-electron integrals according to the basis transformation matrix @param T
     */
    void transform(const Eigen::MatrixXd& T);


/**
 *  Transform the one- and two-electron integrals according to the Jacobi rotation parameters p, q and a given angle theta in radians.
 */
    void rotateJacobi(size_t p, size_t q, double theta);
};


}  // namespace libwint



#endif // LIBWINT_SOBASIS_HPP
