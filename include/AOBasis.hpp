#ifndef LIBWINT_BASIS_HPP
#define LIBWINT_BASIS_HPP


#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "Molecule.hpp"



namespace libwint {


class AOBasis {
private:
    const std::string basisset_name;
    const std::vector<libint2::Atom> atoms;

    // We'd like to keep track if the integrals are calculated already, in order to avoid doing double work
    bool are_calculated_overlap_integrals = false;
    bool are_calculated_nuclear_integrals = false;
    bool are_calculated_kinetic_integrals = false;
    bool are_calculated_electron_repulsion_integrals = false;

    Eigen::MatrixXd S;  // The overlap integrals matrix for the given basis and molecule
    Eigen::MatrixXd V;  // The nuclear integrals matrix for the given basis and molecule
    Eigen::MatrixXd T;  // The kinetic integrals matrix for the given basis and molecule
    Eigen::Tensor<double, 4> g;  // The two-electron repulsion integrals tensor for the given basis and molecule



public:
    // Constructors
    /**
     *  Constructor from a @param: molecule and a @param: basisset_name.
     */
    AOBasis(const libwint::Molecule& molecule, std::string basisset_name);


    // Getters
    Eigen::MatrixXd get_S() const;
    Eigen::MatrixXd get_T() const;
    Eigen::MatrixXd get_V() const;
    Eigen::Tensor<double, 4> get_g() const;


    /**
     *  Calculate and return the number of basis functions in the basis
     */
    size_t calculateNumberOfBasisFunctions() const;

    /**
     *  Calculate and set the overlap integrals, if they haven't been calculated already
     */
    void calculateOverlapIntegrals();

    /**
     *  Calculate and set the kinetic integrals, if they haven't been calculated already
     */
    void calculateKineticIntegrals();

    /**
     *  Calculate and set the nuclear integrals, if they haven't been calculated already
     */
    void calculateNuclearIntegrals();

    /**
     *  Calculate and set the electron repulsion integrals, if they haven't been calculated already
     */
    void calculateElectronRepulsionIntegrals();

    /**
     *  Calculate and set all the integrals, if they haven't been calculated already
     */
    void calculateIntegrals();
};


} // namespace libwint

#endif  // LIBWINT_BASIS_HPP
