#ifndef LIBINT_EIGEN_BASIS_HPP
#define LIBINT_EIGEN_BASIS_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <string>

#include "libint2.hpp"

#include "Molecule.hpp"


namespace Wrapper {

class Basis {
private:
    libint2::BasisSet libint_basis;

public:
    Molecule molecule;
    std::string name;


    // Constructors
    /** Constructor from a molecule and a basis name.
     *
     * @param molecule      Molecule object
     * @param basis_name    string
     */
    Basis(Molecule &molecule, std::string &basis_name);


    // Methods
    size_t nbf();

    Eigen::MatrixXd compute_overlap_integrals();
    Eigen::MatrixXd compute_nuclear_integrals();
    Eigen::MatrixXd compute_kinetic_integrals();
    Eigen::Tensor<double, 4> compute_two_electron_integrals();
};

} // namespace Wrapper

#endif // LIBINT_EIGEN_BASIS_HPP
