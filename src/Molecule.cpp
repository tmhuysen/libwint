#include "Molecule.hpp"


namespace Wrapper {

/** Parses a file name to obtain atoms
 *
 * @param filename
 * @return std::vector<libint2::Atom>
 */
std::vector<libint2::Atom> parse_filename(const std::string& filename) {
    std::cout << "In parse_filename" << std::endl;
    std::cout << "initializing input_file" << std::endl;
    std::ifstream input_file (filename);
    std::cout << "done" << std::endl;
    std::cout << "read_dotxyz" << std::endl;
    return libint2::read_dotxyz (input_file);
}

/** Constructor from a given xyz_filename
 *
 * @param xyz_filename: the path to a .xyz-file that contains the geometry specifications of the molecule.
 *                      IMPORTANT!!! The coordinates of the atoms should be in Angstrom, but LibInt2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
 */
Molecule::Molecule(const std::string& xyz_filename) :
        xyz_filename(xyz_filename)
{
    std::cout << "In Molecule::Molecule" << std::endl;
    this->atoms = parse_filename(this->xyz_filename);
    std::cout << "done parsing filename" << std::endl;
}


unsigned long Molecule::natoms() {
    return this->atoms.size();  // atoms is a std::vector
}

} // namespace Wrapper
