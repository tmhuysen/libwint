#ifndef LIBINT_EIGEN_STUFFEDMOLECULE_HPP
#define LIBINT_EIGEN_STUFFEDMOLECULE_HPP

#include <string>
#include <libint2.hpp>


class StuffedMolecule {
public:

    const std::string xyz_filename;           // Path to a .xyz-file
    std::vector<libint2::Atom> atoms;   // Output of LibInt2's read_dotxyz() function

    /** Constructor from a given xyz_filename
     *
     * @param xyz_filename: the path to a .xyz-file that contains the geometry specifications of the molecule.
     *                      IMPORTANT!!! The coordinates of the atoms should be in Angstrom, but LibInt2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
     */
    StuffedMolecule(const std::string xyz_filename);

    /**
     *
     * @return the number of atoms in the molecule
     */
    unsigned long natoms();

};


#endif //LIBINT_EIGEN_STUFFEDMOLECULE_HPP
