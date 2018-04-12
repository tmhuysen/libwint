#include "SOMullikenBasis.hpp"

namespace libwint {

/*
 *  PRIVATE METHODS
 */

/**
 * Evaluates mulliken operator for two MO's for a given an atomic orbital.
 */
double SOMullikenBasis::evaluateMullikenOperator(size_t molecular_orbital1, size_t molecular_orbital2,
                                                     size_t atomic_orbital) {
    double mulliken_AO_vector_sum= 0;
    for(int i = 0; i<this->K;i++){
        for(int j = 0;j<this->K;j++){
            mulliken_AO_vector_sum += this->S_inverse(atomic_orbital,i)*this->C(j,molecular_orbital2)*this->S(i,j);
        }
    }
    double mulliken_evaluation = 0;
    for(int i = 0; i<this->K;i++){
        mulliken_evaluation += this->C.transpose()(molecular_orbital1,i)*this->S(i,atomic_orbital)*mulliken_AO_vector_sum;

    }
    return mulliken_evaluation;
}


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param atomic orbital instance and a coefficient matrix @param C (i.e. a basis transformation matrix) that links the SO basis to the AO basis
 */

SOMullikenBasis::SOMullikenBasis(const libwint::AOBasis &ao_basis, const Eigen::MatrixXd &C) : SOBasis(ao_basis, C) {
    this->S = ao_basis.get_S();
    this->S_inverse = this->S.inverse();
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

    this->mulliken_matrix = Eigen::MatrixXd::Zero(this->K,this->K);
    for(size_t ao : set_of_AO){
        for(size_t i = 0; i<this->K;i++) {
            double mulliken_evaluation_diagonal = evaluateMullikenOperator(i,i,ao);
            mulliken_matrix(i,i) += mulliken_evaluation_diagonal;
            for (size_t j = 0; j < i; j++) {
                double mulliken_evaluation = evaluateMullikenOperator(i,j,ao)/2 + evaluateMullikenOperator(j,i,ao)/2; // take the hermitian evaluation
                mulliken_matrix(i,j) += mulliken_evaluation;
                mulliken_matrix(j,i) += mulliken_evaluation;
            }
        }
    }
}



/**
 * get_h_SO now returns the one electron value + langrange multiplied corresponding value of the mulliken matrix.
 */
double SOMullikenBasis::get_h_SO(size_t i, size_t j) const {
    return SOBasis::get_h_SO(i, j) + this->lagrange_multiplier*this->mulliken_matrix(i,j);
}


/**
 * Calculate the mulliken population for a set of AO's of a CI a wavefunction (traces the 1RDMs).
 */

double SOMullikenBasis::mullikenPopulationCI(Eigen::MatrixXd &rdm_aa, Eigen::MatrixXd &rdm_bb) {
    double mulliken_population = 0;
    for(size_t i = 0; i<this->K;i++) {
        mulliken_population += mulliken_matrix(i,i)*(rdm_aa(i,i)+rdm_bb(i,i));
    }
    return mulliken_population;
}


double SOMullikenBasis::mullikenPopulationFCI(Eigen::MatrixXd &rdm_aa, Eigen::MatrixXd &rdm_bb) {
    double mulliken_population = 0;
    for(size_t i = 0; i<this->K;i++) {
        for(size_t j = 0; j<this->K;j++){
            mulliken_population += mulliken_matrix(i,j)*(rdm_aa(i,j)+rdm_bb(i,j));
        }
    }
    return mulliken_population;
}

}