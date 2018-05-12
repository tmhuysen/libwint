#ifndef LIBWINT_SOMULLIKENBASIS_HPP
#define LIBWINT_SOMULLIKENBASIS_HPP

#include "SOBasis.hpp"
namespace libwint{
class SOMullikenBasis : public libwint::SOBasis {
private:
    double lagrange_multiplier = 0;
    Eigen::MatrixXd C; //Canonical mat hf
    Eigen::MatrixXd S; //Overlap matrix
    Eigen::MatrixXd mulliken_matrix; //Matrix contains evaluation of the mulliken operator (one electron operator)

    void parseOve(std::string fcidump_filename);
    void parseC(std::string fcidump_filename);

public:
    // Constructors
    /**
     *  Constructor based on a given @param atomic orbital instance and a coefficient matrix @param C (i.e. a basis transformation matrix) that links the SO basis to the AO basis
     */
    SOMullikenBasis(const libwint::AOBasis &ao_basis, const Eigen::MatrixXd &C);

    /**
      *  Constructor based on a given path to an FCIDUMP file
      */
    SOMullikenBasis(std::string fcidump_filename, size_t K);

    //Methods
    /**
     * Calculates the mulliken matrix for a set of AO's, this is the evaluation of the one electron mulliken operator (herm) in MO basis.
     */
    void calculateMullikenMatrix(std::vector<size_t> set_of_AO);

    /**
     * get_h_SO now returns the one electron value + langrange multiplied corresponding value of the mulliken matrix.
     */
    double get_h_SO(size_t i, size_t j) const override;
    Eigen::MatrixXd get_h_SO() const override { return SOBasis::get_h_SO() - this->lagrange_multiplier*this->mulliken_matrix;};
    /**
     * Calculate the mulliken population for a set of AO's of a CI a wavefunction (traces the 1RDMs).
     */
    double mullikenPopulationCI(Eigen::MatrixXd &rdm_aa, Eigen::MatrixXd &rdm_bb);

    void copy(SOMullikenBasis x)  {
        if(this->K != x.K){
            throw std::runtime_error("different base");
        }
        this->h_SO = x.h_SO;
        this->g_SO = x.g_SO;
        this->lagrange_multiplier = x.lagrange_multiplier;
        this->C= x.C;
        this->S= x.S;
        this->mulliken_matrix = x.mulliken_matrix;

    }

    /**
     *  Transform the one- and two-electron integrals (mulliken "integrals" as well) according to the Jacobi rotation parameters p, q and a given angle theta in radians.
     */
    void rotateJacobi(size_t p, size_t q, double theta) override;

    // Setter
    void set_lagrange_multiplier(double lagrange_multiplier) { this->lagrange_multiplier = lagrange_multiplier; }
    // GETTERS
    Eigen::MatrixXd get_mulliken_matrix() { return mulliken_matrix; }
    Eigen::MatrixXd get_C() { return C; }
    Eigen::MatrixXd get_S() { return S; }
    double get_lagrange_multiplier() { return lagrange_multiplier; }





};

}
#endif //LIBWINT_SOMULLIKENBASIS_HPP
