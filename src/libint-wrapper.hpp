//
// Created by Laurent Lemmens on 21/08/17.
//

#ifndef LIBINT_TUTORIAL_LIBINT_WRAPPER_H
#define LIBINT_TUTORIAL_LIBINT_WRAPPER_H

/**
 * Given an operator type, an orbital basis and atoms, calculates the one-body integrals (associated to that operator type)

 * @param opertype: a libint2::Operator (e.g. libint2::Operator::overlap)
 * @param obs:      a libint2::BasisSet object that represents the basis put on the molecule
 * @param atoms:    a std::vector<Atom>
 */
Eigen::MatrixXf compute_1body_integrals(const libint2::Operator& , const libint2::BasisSet& , const std::vector<libint2::Atom> &);


/**
 * Prints the sizes (i.e. the number of basis functions in them) of all shells in a given basis set object.
 *
 * @param obs:  the given basis set object
 */
void print_shell_sizes(const libint2::BasisSet &);

#endif //LIBINT_TUTORIAL_LIBINT_WRAPPER_H
