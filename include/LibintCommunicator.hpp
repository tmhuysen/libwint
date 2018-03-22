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
    static LibintCommunicator& get();


    // Delete the following methods (as required by the singleton class design).
    LibintCommunicator(LibintCommunicator const& libint_communicator) = delete;
    void operator=(LibintCommunicator const& libint_communicator) = delete;


    // Methods
    /**
     *  Calculate the one-body integrals associated to a given @param: operator_type for the given @param: atoms for the basisset with name @param: basisset_name
     */
    Eigen::MatrixXd calculateOneBodyIntegrals(libint2::Operator operator_type, std::string basisset_name, const std::vector<libint2::Atom>& atoms) const;

    /**
     *  Calculate the two-body integrals IN CHEMIST'S NOTATION (11|22) for the given @param: atoms for the basisset with name @param: basisset_name
     */
    Eigen::Tensor<double, 4> calculateTwoBodyIntegrals(std::string basisset_name, const std::vector<libint2::Atom>& atoms) const;
};


} // namespace libwint

#endif // LIBWINT_LIBINTCOMMUNICATOR_HPP
