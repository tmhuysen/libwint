#include "Molecule.hpp"

#include "geometry.hpp"


namespace Wrapper {

/** Parses a file name to obtain atoms
 *
 * @param filename
 * @return std::vector<libint2::Atom>
 */
std::vector<libint2::Atom> parse_filename(const std::string& filename) {
    std::ifstream input_file_stream (filename);
    assert(input_file_stream.good());   // If this assertion fails, we know for sure that we specified a wrong relative path
    return libint2::read_dotxyz (input_file_stream);
}


/** Constructor from a given xyz_filename
 *      The constructed molecule instance corresponds to a neutral atom (i.e. nelec = sum of nucleus charges)
 *
 * @param xyz_filename: the path to a .xyz-file that contains the geometry specifications of the molecule.
 *                      IMPORTANT!!! The coordinates of the atoms should be in Angstrom, but LibInt2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
 */
Molecule::Molecule(const std::string& xyz_filename) :
        xyz_filename(xyz_filename)
{
    this->atoms = parse_filename(this->xyz_filename);

    // For a neutral atom, the number of electrons is equal to the total nucleic charge
    this->nelec = this->nucleic_charge();
}


/** Constructor from a given xyz_filename and a molecular charge
     *      The constructed molecule instance corresponds to an ion:
     *          charge = +1 -> cation (one electron less than the neutral molecule)
     *          charge = 0  -> neutral molecule
     *          charge = -1 -> anion (one electron more than the neutral molecule)
     *
     * @param xyz_filename: the path to a .xyz-file that contains the geometry specifications of the molecule.
     *                      IMPORTANT!!! The coordinates of the atoms should be in Angstrom, but LibInt2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
     */
Molecule::Molecule(const std::string& xyz_filename, int molecular_charge) :
    xyz_filename(xyz_filename)
{
    this->atoms = parse_filename(this->xyz_filename);

    // We're creating an ion here. Since removing an electron increases the charge:
    //  nelec = nucleic_charges - charge
    this->nelec = this->nucleic_charge() - molecular_charge;
}


/** @return the number of atoms in the molecule
 */
size_t Molecule::natoms() {
    return this->atoms.size();  // atoms is a std::vector
}


/** @return the sum of the charges of the nuclei
 */
unsigned Molecule::nucleic_charge() {
    unsigned nucleic_charge = 0;

    for (const auto& atom : this->atoms) {
        nucleic_charge += atom.atomic_number;
    }

    return nucleic_charge;
}


/** @return the internuclear repulsion energy due to the nuclear framework
 *
 */
double Molecule::internuclear_repulsion() {
    double internuclear_repulsion_energy = 0.0;

    auto natoms = this->natoms();
    // Sum over every unique nucleus pair
    for (size_t i = 0; i < natoms; i++) {
        for (size_t j = i + 1; j < natoms; j++ ) {
            const auto atom1 = this->atoms[i];
            const auto atom2 = this->atoms[j];

            // The internuclear repulsion energy (Coulomb) for every nucleus pair is Z1 * Z2 / |R1 - R2|
            internuclear_repulsion_energy += atom1.atomic_number * atom2.atomic_number / Wrapper::distance(atom1, atom2);
        }
    }

    return internuclear_repulsion_energy;
}

} // namespace Wrapper
