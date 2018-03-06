#ifndef LIBWINT_MOLECULE_HPP
#define LIBWINT_MOLECULE_HPP


#include <string>

#include <libint2.hpp>



namespace libwint {


class Molecule {
private:
    const std::vector<libint2::Atom> atoms;
    const size_t N;  // number of electrons


    /**
     *  Parses a @param xyz_filename to @return a std::vector<libint2::Atom>
     */
    std::vector<libint2::Atom> parseXYZFile(std::string xyz_filename) const;



public:
    // Constructors
    /**
     *  Constructor from a given @param xyz_filename
     *      The constructed molecule instance corresponds to a neutral atom (i.e. nelec = sum of nucleus charges)
     *
     *  IMPORTANT!!! The coordinates of the atoms should be input in Angstrom, but libint2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
     */
    Molecule(std::string xyz_filename);

    /**
     *  Constructor from a given @param xyz_filename and a @param molecular_charge
     *      The constructed molecule instance corresponds to an ion:
     *          charge = +1 -> cation (one electron less than the neutral molecule)
     *          charge = 0  -> neutral molecule
     *          charge = -1 -> anion (one electron more than the neutral molecule)
     *
     *  IMPORTANT!!! The coordinates of the atoms should be input in Angstrom, but libint2, which actually processes the .xyz-file, automatically converts to a.u. (bohr).
     */
    Molecule(std::string xyz_filename, int molecular_charge);

    /**
     *  Constructor from a @param atoms: a given std::vector of libint2::Atoms
     *
     *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
     */
    Molecule(const std::vector<libint2::Atom>& atoms);

    /**
     *  Constructor from a @param atoms: a given std::vector of libint2::Atoms and a @param molecular_charge
     *      The constructed molecule instance corresponds to an ion:
     *          charge = +1 -> cation (one electron less than the neutral molecule)
     *          charge = 0  -> neutral molecule
     *          charge = -1 -> anion (one electron more than the neutral molecule)
     *
     *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
     */
    Molecule(const std::vector<libint2::Atom>& atoms, int molecular_charge);


    // Getters
    size_t get_N() const;
    std::vector<libint2::Atom> get_atoms() const;


    // Methods
    /** @return the number of atoms in the molecule
     */
    size_t numberOfAtoms() const;

    /** @return the sum of all the charges of the nuclei
     */
    size_t calculateTotalNucleicCharge() const;

    /** @return the distance between two libint2::Atoms, in Bohr
     */
    double calculateInternuclearDistance(size_t index1, size_t index2) const;

    /** @return the internuclear repulsion energy due to the nuclear framework
     *
     */
    double calculateInternuclearRepulsionEnergy() const;
};


} // namespace libwint

#endif // LIBWINT_MOLECULE_HPP
