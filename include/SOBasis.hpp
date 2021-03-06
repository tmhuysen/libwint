#ifndef LIBWINT_SOBASIS_HPP
#define LIBWINT_SOBASIS_HPP


#include <Eigen/Dense>

#include "AOBasis.hpp"
#include "transformations.hpp"


namespace libwint {


class SOBasis {
protected:
    const size_t K;  // the number of spatial orbitals

    Eigen::MatrixXd h_SO;  // the one-electron integrals (core Hamiltonian) in the spatial orbital basis
    Eigen::Tensor<double, 4> g_SO;  // the two-electron repulsion integrals in the spatial orbital basis


    // Methods
    /**
     *  Parse a given FCIDUMP file for the one- and two-electron integrals
     */
    void parseFCIDUMPFile(std::string fcidump_filename);
    void parseOne(std::string fcidump_filename);
    void parseTwo(std::string fcidump_filename);



public:
    // Constructors
    /**
     *  Constructor based on a given @param atomic orbital instance and a coefficient matrix @param C (i.e. a basis transformation matrix) that links the SO basis to the AO basis
     */
    SOBasis(const libwint::AOBasis& ao_basis, const Eigen::MatrixXd& C);

    explicit SOBasis(size_t K) : K(K){};

    /**
     *  Constructor based on a given path to an FCIDUMP file
     */
    SOBasis(std::string fcidump_filename, size_t K, bool hack = true);

    virtual void copy(SOBasis x) {
        this->h_SO = x.h_SO;
        this->g_SO = x.g_SO;
    };

    // Getters
    const size_t get_K() const { return this->K; }
    virtual Eigen::MatrixXd get_h_SO() const { return this->h_SO; }
    Eigen::Tensor<double, 4> get_g_SO() const { return this->g_SO; }
    virtual double get_h_SO(size_t i, size_t j) const { return this->h_SO(i,j); }
    double get_g_SO(size_t i, size_t j, size_t k, size_t l) const { return this->g_SO(i,j,k,l); }


    // Methods
    /**
     *  Transform the one- and two-electron integrals according to the basis transformation matrix @param T
     */
    void transform(const Eigen::MatrixXd& T);


/**
 *  Transform the one- and two-electron integrals according to the Jacobi rotation parameters p, q and a given angle theta in radians.
 */
    virtual void rotateJacobi(size_t p, size_t q, double theta);
};


}  // namespace libwint



#endif // LIBWINT_SOBASIS_HPP
