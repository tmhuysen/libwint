#include "Molecule.hpp"


namespace Wrapper {

std::vector<libint2::Atom> parse_filename(const std::string& filename) {
    std::ifstream input_file (filename);
    return libint2::read_dotxyz (input_file);
}


Molecule::Molecule(const std::string& xyz_filename) :
        xyz_filename(xyz_filename)
{
    this->atoms = parse_filename(this->xyz_filename);
}


unsigned long Molecule::natoms() {
    return this->atoms.size();  // atoms is a std::vector
}

} // namespace Wrapper
