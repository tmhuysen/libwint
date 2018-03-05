#ifndef LIBWINT_LIBINTCOMMUNICATOR_HPP
#define LIBWINT_LIBINTCOMMUNICATOR_HPP


#include <libint2.hpp>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>



namespace libwint {


/**
 *  A singleton class that takes care of interfacing with the Libint2 (version >2.2.0) C++ API.
 *
 *  Singleton class template from (https://stackoverflow.com/a/1008289).
 */
class LibintCommunicator {
private:
    // Private constructor for the singleton class
    LibintCommunicator();

    // Private destructor: we don't want anyone else to possibly delete the singleton object
    ~LibintCommunicator();



public:
    /**
     *  @return the static singleton instance
     */
    static LibintCommunicator& get() const;


    // Delete the following methods (as required by the singleton class design).
    LibintCommunicator(LibintCommunicator const& libint_communicator) = delete;
    void operator=(LibintCommunicator const& libint_communicator) = delete;


    // Methods
    /**
     * Given an operator type, an orbital basis and atoms, calculates the one-body integrals (associated to that operator type)

     * @param opertype: a libint2::Operator (e.g. libint2::Operator::overlap)
     * @param basis_set:      a libint2::BasisSet object that represents the basis put on the molecule
     * @param atoms:    a std::vector<Atom>

     * @return: an Eigen::MatrixXd storing the integrals
     */
    Eigen::MatrixXd calculateOneBodyIntegrals(libint2::Operator opertype, const libint2::BasisSet& basis_set, const std::vector<libint2::Atom>& atoms) const;

    /**
     * Calculates the two-electron integrals, given an orbital basis and atoms.

     * @param obs:      a libint2::BasisSet object that represents the basis put on the molecule
     * @param atoms:    a std::vector<Atom>

     * @return: an Eigen::Tensor<double, 4> storing the integrals in CHEMIST'S NOTATION (11|22)
     */
    Eigen::Tensor<double, 4> calculateTwoBodyIntegrals(const libint2::BasisSet& obs, const std::vector<libint2::Atom>& atoms) const;
};


} // namespace libwint

#endif // LIBWINT_LIBINTCOMMUNICATOR_HPP
