#include "geometry.hpp"


/** @return the distance between two libint2::Atoms, in Bohr
 */
double libwint::distance(const libint2::Atom& atom1, const libint2::Atom& atom2) {
    return std::sqrt((atom1.x - atom2.x)*(atom1.x - atom2.x) + (atom1.y - atom2.y)*(atom1.y - atom2.y) + (atom1.z - atom2.z)*(atom1.z - atom2.z));
}
