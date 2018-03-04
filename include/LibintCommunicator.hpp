#ifndef LIBWINT_LIBINTCOMMUNICATOR_HPP
#define LIBWINT_LIBINTCOMMUNICATOR_HPP


#include <libwint.hpp>
#include <Eigen/Dense>



namespace libwint {


/**
 *  A singleton class that takes care of interfacing with the Libint2 (version >2.2.0) C++ API.
 *
 *  Singleton class template from (https://stackoverflow.com/a/1008289).
 */
class LibintCommunicator {
private:
    // Constructor
    LibintCommunicator();

    // Destructor: we don't want anyone else to possibly delete the singleton object
    ~LibintCommunicator();


    /**
     * Given an operator type, an orbital basis and atoms, calculates the one-body integrals (associated to that operator type)

     * @param opertype: a libint2::Operator (e.g. libint2::Operator::overlap)
     * @param obs:      a libint2::BasisSet object that represents the basis put on the molecule
     * @param atoms:    a std::vector<Atom>

     * @return: an Eigen::MatrixXd storing the integrals
     */
    Eigen::MatrixXd computeOneBodyIntegrals(const libint2::Operator& opertype, const libint2::BasisSet& obs, const std::vector<libint2::Atom>& atoms);

    /**
     * Calculates the two-electron integrals, given an orbital basis and atoms.

     * @param obs:      a libint2::BasisSet object that represents the basis put on the molecule
     * @param atoms:    a std::vector<Atom>

     * @return: an Eigen::Tensor<double, 4> storing the integrals in CHEMIST'S NOTATION (11|22)
     */
    Eigen::Tensor<double, 4> computeTwoBodyIntegrals(const libint2::BasisSet& obs, const std::vector<libint2::Atom>& atoms);



public:
    /**
     *  @return the static singleton instance
     */
    static LibintCommunicator& get();


    // Delete the following methods (as required by the singleton class design).
    LibintCommunicator(LibintCommunicator const& libint_communicator) = delete;
    void operator=(LibintCommunicator const& libint_communicator) = delete;
};


} // namespace libwint

#endif // LIBWINT_LIBINTCOMMUNICATOR_HPP
