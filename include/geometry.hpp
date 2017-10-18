#ifndef LIBINT_WRAPPER_GEOMETRY_HPP
#define LIBINT_WRAPPER_GEOMETRY_HPP

#include <libint2.hpp>


namespace Wrapper{

/** @return the distance between two libint2::Atoms
 */
double distance(libint2::Atom& atom1, libint2::Atom& atom2);


}


#endif // LIBINT_WRAPPER_GEOMETRY_HPP
