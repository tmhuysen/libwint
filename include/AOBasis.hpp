#ifndef LIBWINT_BASIS_HPP
#define LIBWINT_BASIS_HPP


#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "Molecule.hpp"



namespace libwint {


class AOBasis {
private:
    // We'd like to keep track if the integrals are calculated already, in order to avoid doing double work
    bool are_calculated_overlap_integrals = false;
    bool are_calculated_nuclear_integrals = false;
    bool are_calculated_kinetic_integrals = false;
    bool are_calculated_electron_repulsion_integrals = false;

    const Molecule& molecule;
    const std::string name;
    const libint2::BasisSet& libint_basis;

    Eigen::MatrixXd S;  // The overlap integrals matrix for the given basis and molecule
    Eigen::MatrixXd V;  // The nuclear integrals matrix for the given basis and molecule
    Eigen::MatrixXd T;  // The kinetic integrals matrix for the given basis and molecule
    Eigen::Tensor<double, 4> g;  // The two-electron repulsion integrals tensor for the given basis and molecule



public:
    // Constructors
    /** Constructor from a molecule and a basis name.
     *
     * @param molecule      Molecule object
     * @param basis_name    string
     */
    AOBasis(Molecule& molecule, std::string basis_name);


    // Getters
    libwint::Molecule get_molecule() const;
    Eigen::MatrixXd get_S() const;
    Eigen::MatrixXd get_T() const;
    Eigen::MatrixXd get_V() const;
    Eigen::Tensor<double, 4> get_g() const;


    /** Calculate and return the number of basis functions in the basis
     */
    size_t calculateNumberOfBasisFunctions() const;

    /** Calculate and set the overlap integrals
     */
    void calculateOverlapIntegrals();

    /** Calculate and set the kinetic integrals
    */
    void calculateKineticIntegrals();

    /** Calculate and set the nuclear integrals
    */
    void calculateNuclearIntegrals();

    /** Calculate and set the electron repulsion integrals
    */
    void calculateElectronRepulsionIntegrals();

    /** Calculate and set all the integrals
     */
    void calculateIntegrals();
};


} // namespace libwint

#endif // LIBWINT_BASIS_HPP
