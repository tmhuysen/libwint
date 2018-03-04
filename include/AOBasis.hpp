#ifndef LIBWINT_BASIS_HPP
#define LIBWINT_BASIS_HPP


#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "LibintCommunicator.hpp
#include "Molecule.hpp"



namespace libwint {


class AOBasis {
private:
    // We'd like to keep track if the integrals are calculated already, in order to avoid doing double work
    bool are_computed_overlap_integrals = false;
    bool are_computed_nuclear_integrals = false;
    bool are_computed_kinetic_integrals = false;
    bool are_computed_tei = false;

    libint2::BasisSet& libint_basis;
    const std::string name;
    Molecule& molecule;

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
    AOBasis(Molecule& molecule, const std::string& basis_name);


    /** Calculate and return the number of basis functions in the basis
     */
    size_t nbf();

    /** Calculate and set the overlap integrals
     */
    void compute_overlap_integrals();

    /** Calculate and set the kinetic integrals
    */
    void compute_kinetic_integrals();

    /** Calculate and set the nuclear integrals
    */
    void compute_nuclear_integrals();

    /** Calculate and set the kinetic integrals
    */
    void compute_two_electron_integrals();

    /** Calculate and set all the integrals
     */
    void compute_integrals();
};

} // namespace libwint

#endif // LIBWINT_BASIS_HPP
