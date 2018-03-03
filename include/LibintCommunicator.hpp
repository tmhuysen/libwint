#ifndef LIBWINT_LIBINTCOMMUNICATOR_HPP
#define LIBWINT_LIBINTCOMMUNICATOR_HPP

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
