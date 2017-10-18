#include "geometry.hpp"


namespace Wrapper {

/** @return the distance between two libint2::Atoms
 */
double distance(libint2::Atom& atom1, libint2::Atom& atom2) {
    return std::sqrt((atom1.x - atom2.x)*(atom1.x - atom2.x) + (atom1.y - atom2.y)*(atom1.y - atom2.y) + (atom1.z - atom2.z)*(atom1.z - atom2.z));
}


}