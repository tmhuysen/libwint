#include "SOMullikenBasis.hpp"

namespace libwint {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param atomic orbital instance and a coefficient matrix @param C (i.e. a basis transformation matrix) that links the SO basis to the AO basis
 */

SOMullikenBasis::SOMullikenBasis(const libwint::AOBasis &ao_basis, const Eigen::MatrixXd &C) : SOBasis(ao_basis, C) {
    this->S = ao_basis.get_S();
    this->C = C;
    this->mulliken_matrix = Eigen::MatrixXd::Zero(this->K, this->K);
}


/*
 * Public methods
 */

/**
 * Calculates the mulliken matrix for a set of AO's, this is the evaluation of the one electron mulliken (herm) in MO basis.
 */

void SOMullikenBasis::calculateMullikenMatrix(std::vector<size_t> set_of_AO) {

    Eigen::MatrixXd p_a = Eigen::MatrixXd::Zero(this->K,this->K);
    Eigen::MatrixXd Ct = Eigen::MatrixXd(this->C.transpose());
    for(size_t ao : set_of_AO) { p_a(ao, ao) = 1; }
    this->mulliken_matrix = (Ct*p_a*this->S*this->C + Ct*this->S*p_a*this->C)/2;





}


/**
 * get_h_SO now returns the one electron value + langrange multiplied corresponding value of the mulliken matrix.
 */
double SOMullikenBasis::get_h_SO(size_t i, size_t j) const {
    return SOBasis::get_h_SO(i, j) - this->lagrange_multiplier*this->mulliken_matrix(i,j);
}


/**
 * Calculate the mulliken population for a set of AO's of a CI a wavefunction (traces the 1RDMs).
 */

double SOMullikenBasis::mullikenPopulationCI(Eigen::MatrixXd &rdm_aa, Eigen::MatrixXd &rdm_bb) {
    double mulliken_population =(mulliken_matrix*(rdm_aa+rdm_bb)).trace();
    return mulliken_population;
}

}  // namespace libwint