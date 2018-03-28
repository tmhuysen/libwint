#ifndef LIBWINT_SOMULLIKENBASIS_HPP
#define LIBWINT_SOMULLIKENBASIS_HPP

#include "SOBasis.hpp"
namespace libwint{
class SOMullikenBasis : public libwint::SOBasis {
private:
    double lagrange_multiplier = 0;

    Eigen::MatrixXd C; //Canonical mat hf
    Eigen::MatrixXd S; //Overlap matrix
    Eigen::MatrixXd S_inverse; //Overlap matrix
    Eigen::MatrixXd mulliken_matrix; //Matrix contains evaluation of the mulliken operator (one electron operator)


    // Methods
    /**
     * Evaluates mulliken operator for two MO's for a given an atomic orbital.
     */
    double evaluateMullikenOperator(size_t molecular_orbital1, size_t molecular_orbital2, size_t atomic_orbital);


public:
    // Constructors
    /**
     *  Constructor based on a given @param atomic orbital instance and a coefficient matrix @param C (i.e. a basis transformation matrix) that links the SO basis to the AO basis
     */
    SOMullikenBasis(const libwint::AOBasis &ao_basis, const Eigen::MatrixXd &C);


    //Methods
    /**
     * Calculates the mulliken matrix for a set of AO's, this is the evaluation of the one electron mulliken (herm) in MO basis.
     */
    void calculateMullikenMatrix(std::vector<size_t> set_of_AO);

    /**
     * get_h_SO now returns the one electron value + langrange multiplied corresponding value of the mulliken matrix.
     */
    double get_h_SO(size_t i, size_t j) const override;

    /**
     * Calculate the mulliken population for a set of AO's of a CI a wavefunction (traces the 1RDMs).
     */
    double mullikenPopulationCI(Eigen::MatrixXd &rdm_aa, Eigen::MatrixXd &rdm_bb);


    // Setters
    void set_lagrange_multiplier(double lagrange_multiplier) { this->lagrange_multiplier = lagrange_multiplier; }
};

}
#endif //LIBWINT_SOMULLIKENBASIS_HPP
