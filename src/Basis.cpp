#include "Basis.hpp"

#include "integrals.hpp"


/** Constructor from a molecule and a basis name.
 *
 * @param molecule      Molecule object
 * @param basis_name    string
 */
libwrp::Basis::Basis(Molecule& molecule, const std::string& basis_name) :
        molecule(molecule), name(basis_name) {
    // Constructing the basis also constructs the associated libint2::BasisSet object
    libint2::BasisSet libint_basis(this->name, this->molecule.atoms);
    this->libint_basis = libint_basis;
}


/** Calculate and return the number of basis functions in the basis
 */
size_t libwrp::Basis::nbf() {
    return static_cast<size_t>(this->libint_basis.nbf());
}


/** Calculate and set the overlap integrals
 *
 *      If the overlap integrals have already been calculated, print an error message
 */
void libwrp::Basis::compute_overlap_integrals() {

    if (!this->are_computed_overlap_integrals) {
        this->S = compute_1body_integrals(libint2::Operator::overlap, this->libint_basis, this->molecule.atoms);

        this->are_computed_overlap_integrals = true;
    } else {
        std::cout << "The overlap integrals have already been calculated in this basis ..." << std::endl;
    }

}


/** Calculate and set the kinetic integrals
*/
void libwrp::Basis::compute_kinetic_integrals() {

    if (!this->are_computed_kinetic_integrals) {
        this -> T = compute_1body_integrals(libint2::Operator::kinetic, this->libint_basis, this->molecule.atoms);

        this->are_computed_kinetic_integrals = true;
    } else {
        std::cout << "The kinetic integrals have already been calculated in this basis ..." << std::endl;
    }

}


/** Calculate and set the nuclear integrals
*/
void libwrp::Basis::compute_nuclear_integrals() {

    if (!this->are_computed_nuclear_integrals) {
        this-> V = compute_1body_integrals(libint2::Operator::nuclear, this->libint_basis, this->molecule.atoms);

        this->are_computed_nuclear_integrals = true;
    } else {
        std::cout << "The nuclear integrals have already been calculated in this basis ..." << std::endl;
    }

}


/** Calculate and set the kinetic integrals
*/
void libwrp::Basis::compute_two_electron_integrals() {

    if (!this->are_computed_tei) {
        this -> tei = compute_2body_integrals(this->libint_basis, this->molecule.atoms);

        this->are_computed_tei = true;
    } else {
        std::cout << "The two-electron integrals have already been calculated in this basis ..." << std::endl;
    }

};


/** Calculate and set all the integrals
 */
void libwrp::Basis::compute_integrals() {
    compute_overlap_integrals();
    compute_kinetic_integrals();
    compute_nuclear_integrals();
    compute_two_electron_integrals();
}
