#ifndef LIBWINT_GEOMETRY_HPP
#define LIBWINT_GEOMETRY_HPP

#include <libint2.hpp>


namespace libwint {

/** @return the distance between two libint2::Atoms, in Bohr
 */
double distance(const libint2::Atom& atom1, const libint2::Atom& atom2);


} // namespace libwint


#endif // LIBWINT_GEOMETRY_HPP
