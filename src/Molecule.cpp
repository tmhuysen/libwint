#include "Molecule.hpp"

#include <boost/filesystem.hpp>



namespace libwint {


/*
 *  PRIVATE METHODS
 */

/**
 *  Parses a @param filename to @return a std::vector<libint2::Atom>
 */
std::vector<libint2::Atom> Molecule::parseXYZFilename(const std::string filename) const {

    // If the filename doesn't end with .xyz, we assume that the user supplied a wrong file
    boost::filesystem::path path (filename);
    std::string extension = path.filename().string();
    if (!(extension == ".xyz")) {
        throw std::runtime_error("You did not provide a .xyz file name");
    }

    // If the filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
    std::ifstream input_file_stream (filename);
    if (!input_file_stream.good()) {
        throw std::runtime_error("The provided .xyz file name is illegible. Maybe you specified a wrong path?");
    }

    return libint2::read_dotxyz (input_file_stream);  // can't make a reference because that's how libint2 is implemented
}



/*
 *  CONSTRUCTORS
 */

/** Constructor from a given xyz_filename
 *      The constructed molecule instance corresponds to a neutral atom (i.e. nelec = sum of nucleus charges)
 *
 * @param xyz_filename: the path to a .xyz-file that contains the geometry specifications of the molecule.
 *                      IMPORTANT!!! The coordinates of the atoms should be in Angstrom, but LibInt2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
 */
Molecule::Molecule(const std::string& xyz_filename) :
        atoms (this->parseXYZFilename(xyz_filename)),
        N (this->calculateTotalNucleicCharge())
{}


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
        atoms (this->parseXYZFilename(xyz_filename)),
        N (this->calculateTotalNucleicCharge() - molecular_charge)  // we're possibly creating an ion, so N = charges of nuclei - total molecular charge
{}


/**
 *  Constructor from a @param atoms: a given std::vector of libint2::Atoms
 *
 *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
 */
Molecule::Molecule(const std::vector<libint2::Atom>& atoms) :
    atoms (atoms),
    N (this->calculateTotalNucleicCharge())
{}


/**
 *  Constructor from a @param atoms: a given std::vector of libint2::Atoms and a @param molecular_charge
 *      The constructed molecule instance corresponds to an ion:
 *          charge = +1 -> cation (one electron less than the neutral molecule)
 *          charge = 0  -> neutral molecule
 *          charge = -1 -> anion (one electron more than the neutral molecule)
 *
 *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
 */
Molecule::Molecule(const std::vector<libint2::Atom>& atoms, int molecular_charge) :
    atoms (atoms),
    N (this->calculateTotalNucleicCharge() - molecular_charge)
{}



/*
 *  GETTERS
 */

size_t Molecule::get_N() const { return this->N; }



/*
 *  PUBLIC METHODS
 */

/** @return the number of atoms in the molecule
 */
size_t Molecule::numberOfAtoms() const { return this->atoms.size(); }  // atoms is a std::vector


/** @return the sum of all the charges of the nuclei
 */
size_t Molecule::calculateTotalNucleicCharge() const {
    size_t nucleic_charge = 0;

    for (const auto& atom : this->atoms) {
        nucleic_charge += atom.atomic_number;
    }

    return nucleic_charge;
}


/** @return the distance between two libint2::Atoms, in Bohr
 */
double Molecule::calculateInternuclearDistance(size_t index1, size_t index2) const {

    const libint2::Atom atom1 = this->atoms[index1];
    const libint2::Atom atom2 = this->atoms[index2];

    return std::sqrt((atom1.x - atom2.x)*(atom1.x - atom2.x) + (atom1.y - atom2.y)*(atom1.y - atom2.y) + (atom1.z - atom2.z)*(atom1.z - atom2.z));
}


/** @return the internuclear repulsion energy due to the nuclear framework
 *
 */
double Molecule::calculateInternuclearRepulsionEnergy() const {
    double internuclear_repulsion_energy = 0.0;

    auto natoms = this->numberOfAtoms();
    // Sum over every unique nucleus pair
    for (size_t i = 0; i < natoms; i++) {
        for (size_t j = i + 1; j < natoms; j++ ) {
            const auto atom1 = this->atoms[i];
            const auto atom2 = this->atoms[j];

            // The internuclear repulsion energy (Coulomb) for every nucleus pair is Z1 * Z2 / |R1 - R2|
            internuclear_repulsion_energy += atom1.atomic_number * atom2.atomic_number / this->calculateInternuclearDistance(i, j);
        }
    }

    return internuclear_repulsion_energy;
}


}  // namespace libwint
