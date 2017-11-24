#include "Basis.hpp"

#include "integrals.hpp"


/** Constructor from a molecule and a basis name.
 *
 * @param molecule      Molecule object
 * @param basis_name    string
 */
Wrapper::Basis::Basis(Molecule& molecule, const std::string& basis_name) :
        molecule(molecule), name(basis_name) {
    // Constructing the basis also constructs the associated libint2::BasisSet object
    libint2::BasisSet libint_basis(this->name, this->molecule.atoms);
    this->libint_basis = libint_basis;
}


/** Calculate and return the number of basis functions in the basis
 */
size_t Wrapper::Basis::nbf() {
    return static_cast<size_t>(this->libint_basis.nbf());
}


/** Calculate and set the overlap integrals
 */
void Wrapper::Basis::compute_overlap_integrals() {
    this->S = compute_1body_integrals(libint2::Operator::overlap, this->libint_basis, this->molecule.atoms);
}


/** Calculate and set the kinetic integrals
*/
void Wrapper::Basis::compute_kinetic_integrals() {
    this -> T = compute_1body_integrals(libint2::Operator::kinetic, this->libint_basis, this->molecule.atoms);
}


/** Calculate and set the nuclear integrals
*/
void Wrapper::Basis::compute_nuclear_integrals() {
    this-> V = compute_1body_integrals(libint2::Operator::nuclear, this->libint_basis, this->molecule.atoms);
}


/** Calculate and set the kinetic integrals
*/
void Wrapper::Basis::compute_two_electron_integrals() {
    this -> tei = compute_2body_integrals(this->libint_basis, this->molecule.atoms);
};
