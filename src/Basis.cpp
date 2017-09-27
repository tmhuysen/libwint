#include "Basis.hpp"

#include "integrals.hpp"


namespace Wrapper {

/** Constructor from a molecule and a basis name.
 *
 * @param molecule      Molecule object
 * @param basis_name    string
 */
Basis::Basis(Molecule &molecule, std::string &basis_name) :
        molecule(molecule), name(basis_name) {
    // Constructing the basis also constructs the associated libint2::BasisSet object
    libint2::BasisSet libint_basis(this->name, this->molecule.atoms);
    this->libint_basis = libint_basis;
}

size_t Basis::nbf() {
    return this->libint_basis.nbf();
}

Eigen::MatrixXd Basis::compute_overlap_integrals() {
    return compute_1body_integrals(libint2::Operator::overlap, this->libint_basis, this->molecule.atoms);
}

Eigen::MatrixXd Basis::compute_kinetic_integrals() {
    return compute_1body_integrals(libint2::Operator::kinetic, this->libint_basis, this->molecule.atoms);
}

Eigen::MatrixXd Basis::compute_nuclear_integrals() {
    return compute_1body_integrals(libint2::Operator::nuclear, this->libint_basis, this->molecule.atoms);
}

Eigen::Tensor<double, 4> Basis::compute_two_electron_integrals() {
    return compute_2body_integrals(this->libint_basis, this->molecule.atoms);
};

} // namespace Wrapper
