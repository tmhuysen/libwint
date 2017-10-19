#ifndef LIBWRP_BASIS_HPP
#define LIBWRP_BASIS_HPP

#include "Molecule.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


namespace Wrapper {

class Basis {
private:
    libint2::BasisSet libint_basis;

public:
    Molecule molecule;
    const std::string name;


    // Constructors
    /** Constructor from a molecule and a basis name.
     *
     * @param molecule      Molecule object
     * @param basis_name    string
     */
    Basis(Molecule& molecule, const std::string& basis_name);


    // Methods
    size_t nbf();

    Eigen::MatrixXd compute_overlap_integrals();
    Eigen::MatrixXd compute_nuclear_integrals();
    Eigen::MatrixXd compute_kinetic_integrals();
    Eigen::Tensor<double, 4> compute_two_electron_integrals();
};

} // namespace Wrapper

#endif // LIBWRP_BASIS_HPP
