#ifndef LIBINT_WRAPPER_MOLECULE_HPP
#define LIBINT_WRAPPER_MOLECULE_HPP

#include <libint2.hpp>
#include <string>


namespace Wrapper {
/** Parses a file name to obtain atoms
 *
 * @param filename
 * @return std::vector<libint2::Atom>
 */
std::vector<libint2::Atom> parse_filename(const std::string &filename);


class Molecule {
public:
    const std::string xyz_filename;     // Path to a .xyz-file
    std::vector<libint2::Atom> atoms;   // Output of LibInt2's read_dotxyz() function


    // Constructors
    /** Constructor from a given xyz_filename
     *
     * @param xyz_filename: the path to a .xyz-file that contains the geometry specifications of the molecule.
     *                      IMPORTANT!!! The coordinates of the atoms should be in Angstrom, but LibInt2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
     */
    Molecule(const std::string &xyz_filename);


    // Methods
    /**
     *
     * @return the number of atoms in the molecule
     */
    unsigned long natoms();
};

} // namespace Wrapper

#endif // LIBINT_WRAPPER_MOLECULE_HPP
