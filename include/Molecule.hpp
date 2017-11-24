#ifndef LIBWRP_MOLECULE_HPP
#define LIBWRP_MOLECULE_HPP

#include <libint2.hpp>
#include <string>


namespace libwrp {
/** Parses a file name to obtain atoms
 *
 * @param filename
 * @return std::vector<libint2::Atom>
 */
std::vector<libint2::Atom> parse_filename(const std::string& filename);


class Molecule {
public:
    const std::string xyz_filename;     // Path to a .xyz-file
    std::vector<libint2::Atom> atoms;   // A std::vector of libint2::Atoms
                                        //      A libint2::Atom is just a struct with data fields charge, x, y ,z
                                        //      Note that x, y, and z are all in Bohr, but the .xyz-file should specify them in Angstrom

    unsigned nelec;                     // The number of electrons in the molecule


    // Constructors
    /** Constructor from a given xyz_filename
     *      The constructed molecule instance corresponds to a neutral atom (i.e. nelec = sum of nucleus charges)
     *
     * @param xyz_filename: the path to a .xyz-file that contains the geometry specifications of the molecule.
     *                      IMPORTANT!!! The coordinates of the atoms should be in Angstrom, but LibInt2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
     */
    Molecule(const std::string& xyz_filename);

    /** Constructor from a given xyz_filename and a molecular charge
     *      The constructed molecule instance corresponds to an ion:
     *          charge = +1 -> cation (one electron less than the neutral molecule)
     *          charge = 0  -> neutral molecule
     *          charge = -1 -> anion (one electron more than the neutral molecule)
     *
     * @param xyz_filename: the path to a .xyz-file that contains the geometry specifications of the molecule.
     *                      IMPORTANT!!! The coordinates of the atoms should be in Angstrom, but LibInt2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
     */
    Molecule(const std::string& xyz_filename, int molecular_charge);


    // Methods
    /** @return the number of atoms in the molecule
     */
    size_t natoms();

    /** @return the sum of the charges of the nuclei
     */
    unsigned nucleic_charge();

    /** @return the internuclear repulsion energy due to the nuclear framework
     *
     */
    double internuclear_repulsion();
};

} // namespace libwrp

#endif // LIBWRP_MOLECULE_HPP
