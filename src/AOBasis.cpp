#include "AOBasis.hpp"

#include "LibintCommunicator.hpp"



namespace libwint {

/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor from a @param: molecule and a @param: basis name.
 */
AOBasis::AOBasis(const Molecule& molecule, std::string basisset_name) :
        basisset_name (basisset_name),
        atoms (molecule.get_atoms())
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

/**
 * Calculate and return the number of basis functions in the basis
 */
size_t AOBasis::calculateNumberOfBasisFunctions() const {

    // Find out if any of the integrals have been calculated yet and return the dimension
    if (this->are_calculated_overlap_integrals) {
        return static_cast<size_t>(this->S.cols());
    }
    else if (this->are_calculated_kinetic_integrals) {
        return static_cast<size_t>(this->T.cols());
    }
    else if (this->are_calculated_nuclear_integrals) {
        return static_cast<size_t>(this->V.cols());
    }
    else if (this->are_calculated_electron_repulsion_integrals) {
        return static_cast<size_t>(this->g.dimension(0));
    }
    else {
        throw std::runtime_error("None of the integrals have been calculated yet, so I can't figure out how many basis functions there are.");
    }
}


/**
 *  Calculate and set the overlap integrals, if they haven't been calculated already
 */
void AOBasis::calculateOverlapIntegrals() {

    if (!this->are_calculated_overlap_integrals) {
        this->S = libwint::LibintCommunicator::get().calculateOneBodyIntegrals(libint2::Operator::overlap, this->basisset_name, this->atoms);
        this->are_calculated_overlap_integrals = true;
    } else {
        std::cout << "The overlap integrals have already been calculated in this basis ..." << std::endl;
    }
}


/**
 *  Calculate and set the kinetic integrals, if they haven't been calculated already
 */
void AOBasis::calculateKineticIntegrals() {

    if (!this->are_calculated_kinetic_integrals) {
        this -> T = libwint::LibintCommunicator::get().calculateOneBodyIntegrals(libint2::Operator::kinetic, this->basisset_name, this->atoms);
        this->are_calculated_kinetic_integrals = true;
    } else {
        std::cout << "The kinetic integrals have already been calculated in this basis ..." << std::endl;
    }
}


/**
 *  Calculate and set the nuclear integrals, if they haven't been calculated already
 */
void AOBasis::calculateNuclearIntegrals() {

    if (!this->are_calculated_nuclear_integrals) {
        this-> V = libwint::LibintCommunicator::get().calculateOneBodyIntegrals(libint2::Operator::nuclear, this->basisset_name, this->atoms);
        this->are_calculated_nuclear_integrals = true;
    } else {
        std::cout << "The nuclear integrals have already been calculated in this basis ..." << std::endl;
    }
}


/**
 *  Calculate and set the kinetic integrals, if they haven't been calculated already
 */
void AOBasis::calculateElectronRepulsionIntegrals() {

    if (!this->are_calculated_electron_repulsion_integrals) {
        this -> g = libwint::LibintCommunicator::get().calculateTwoBodyIntegrals(this->basisset_name, this->atoms);
        this->are_calculated_electron_repulsion_integrals = true;
    } else {
        std::cout << "The two-electron repulsion integrals have already been calculated in this basis ..." << std::endl;
    }
};


/**
 *  Calculate and set all the integrals, if they haven't been calculated already
 */
void AOBasis::calculateIntegrals() {
    calculateOverlapIntegrals();
    calculateKineticIntegrals();
    calculateNuclearIntegrals();
    calculateElectronRepulsionIntegrals();
}


}  // namespace libwint
