#ifndef LIBWRP_BASIS_HPP
#define LIBWRP_BASIS_HPP

#include "Molecule.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


namespace libwrp {

class Basis {
private:
    libint2::BasisSet libint_basis;

    bool are_computed_overlap_integrals = false;
    bool are_computed_nuclear_integrals = false;
    bool are_computed_kinetic_integrals = false;
    bool are_computed_tei = false;

public:
    Molecule& molecule;         // Make a reference to a Molecule object
    const std::string name;

    Eigen::MatrixXd S;  // The overlap integrals matrix for the given basis and molecule
    Eigen::MatrixXd V;  // The nuclear integrals matrix for the given basis and molecule
    Eigen::MatrixXd T;  // The kinetic integrals matrix for the given basis and molecule
    Eigen::Tensor<double, 4> tei;  // The two-electron repulsion integrals tensor for the given basis and molecule


    // Constructors
    /** Constructor from a molecule and a basis name.
     *
     * @param molecule      Molecule object
     * @param basis_name    string
     */
    Basis(Molecule& molecule, const std::string& basis_name);


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

} // namespace libwrp

#endif // LIBWRP_BASIS_HPP
