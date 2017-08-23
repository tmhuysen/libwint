//
// Created by Laurent Lemmens on 22/08/17.
//

#include "StuffedMolecule.hpp"


Molecule::Molecule(const std::string xyz_filename) : xyz_filename(xyz_filename) {
    // FIXME: Raise error if doesn't end with .xyz
    std::ifstream input_file (this->xyz_filename);
    this->atoms = libint2::read_dotxyz(input_file);
}


unsigned long Molecule::natoms() {
    return this->atoms.size();  // atoms is a std::vector
}