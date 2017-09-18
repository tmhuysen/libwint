//
// Created by Laurent Lemmens on 21/08/17.
//

#ifndef LIBINT_TUTORIAL_LIBINT_WRAPPER_H
#define LIBINT_TUTORIAL_LIBINT_WRAPPER_H

#include <libint2.hpp>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>


/**
 * Given an operator type, an orbital basis and atoms, calculates the one-body integrals (associated to that operator type).
 */
Eigen::MatrixXd compute_1body_integrals(const libint2::Operator&, const libint2::BasisSet&, const std::vector<libint2::Atom>&);

/**
 * Calculates the two-electron integrals, given an orbital basis and atoms.
 */
Eigen::Tensor<double, 4> compute_2body_integrals(const libint2::BasisSet&, const std::vector<libint2::Atom>&);


/**
 * Prints the sizes (i.e. the number of basis functions in them) of all shells in a given basis set object.
 */
void print_shell_sizes(const libint2::BasisSet &);

#endif //LIBINT_TUTORIAL_LIBINT_WRAPPER_H
