#include "AOBasis.hpp"

#include "LibintCommunicator.hpp"



namespace libwint {

/*
 *  CONSTRUCTORS
 */

/** Constructor from a molecule and a basis name.
 *
 * @param molecule      Molecule object
 * @param basis_name    string
 */
AOBasis::AOBasis(Molecule& molecule, const std::string basis_name) :
        molecule(molecule),
        name(basis_name),
        libint_basis(libint2::BasisSet (this->name, this->molecule.get_atoms()))
{}



/*
 *  GETTERS
 */

Eigen::MatrixXd AOBasis::get_S() const {

    if (!this->are_calculated_overlap_integrals) {
        throw std::logic_error("You haven't calculated the overlap integrals yet and are trying to access them.");
    } else {
        return this->S;
    }
}

Eigen::MatrixXd AOBasis::get_T() const {

    if (!this->are_calculated_kinetic_integrals) {
        throw std::logic_error("You haven't calculated the kinetic integrals yet and are trying to access them.");
    } else {
        return this->T;
    }
}

Eigen::MatrixXd AOBasis::get_V() const {

    if (!this->are_calculated_nuclear_integrals) {
        throw std::logic_error("You haven't calculated the nuclear integrals yet and are trying to access them.");
    } else {
        return this->V;
    }
}

Eigen::Tensor<double, 4> AOBasis::get_g() const {

    if (!this->are_calculated_electron_repulsion_integrals) {
        throw std::logic_error("You haven't calculated the electron repulsion integrals yet and are trying to access them.");
    } else {
        return this->g;
    }
};



/*
 *  PUBLIC METHODS
 */

/** Calculate and return the number of basis functions in the basis
 */
size_t AOBasis::calculateNumberOfBasisFunctions() const {
    return static_cast<size_t>(this->libint_basis.nbf());
}


/** Calculate and set the overlap integrals
 *
 *      If the overlap integrals have already been calculated, print an error message
 */
void AOBasis::calculateOverlapIntegrals() {

    if (!this->are_calculated_overlap_integrals) {
        this->S = libwint::LibintCommunicator::get().calculateOneBodyIntegrals(libint2::Operator::overlap, this->libint_basis, this->molecule.get_atoms());
        this->are_calculated_overlap_integrals = true;
    } else {
        std::cout << "The overlap integrals have already been calculated in this basis ..." << std::endl;
    }

}


/** Calculate and set the kinetic integrals
*/
void AOBasis::calculateKineticIntegrals() {

    if (!this->are_calculated_kinetic_integrals) {
        this -> T = libwint::LibintCommunicator::get().calculateOneBodyIntegrals(libint2::Operator::kinetic, this->libint_basis, this->molecule.get_atoms());
        this->are_calculated_kinetic_integrals = true;
    } else {
        std::cout << "The kinetic integrals have already been calculated in this basis ..." << std::endl;
    }

}


/** Calculate and set the nuclear integrals
*/
void AOBasis::calculateNuclearIntegrals() {

    if (!this->are_calculated_nuclear_integrals) {
        this-> V = libwint::LibintCommunicator::get().calculateOneBodyIntegrals(libint2::Operator::nuclear, this->libint_basis, this->molecule.get_atoms());
        this->are_calculated_nuclear_integrals = true;
    } else {
        std::cout << "The nuclear integrals have already been calculated in this basis ..." << std::endl;
    }

}


/** Calculate and set the kinetic integrals
*/
void AOBasis::calculateElectronRepulsionIntegrals() {

    if (!this->are_calculated_electron_repulsion_integrals) {
        this -> g = libwint::LibintCommunicator::get().calculateTwoBodyIntegrals(this->libint_basis, this->molecule.get_atoms());
        this->are_calculated_electron_repulsion_integrals = true;
    } else {
        std::cout << "The two-electron integrals have already been calculated in this basis ..." << std::endl;
    }

};


/** Calculate and set all the integrals
 */
void AOBasis::calculateIntegrals() {
    calculateOverlapIntegrals();
    calculateKineticIntegrals();
    calculateNuclearIntegrals();
    calculateElectronRepulsionIntegrals();
}


}  // namespace libwint
