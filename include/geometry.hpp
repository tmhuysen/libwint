#ifndef LIBWRP_GEOMETRY_HPP
#define LIBWRP_GEOMETRY_HPP

#include <libint2.hpp>


namespace libwrp {

/** @return the distance between two libint2::Atoms, in Bohr
 */
double distance(const libint2::Atom& atom1, const libint2::Atom& atom2);


} // namespace libwrp


#endif // LIBWRP_GEOMETRY_HPP
