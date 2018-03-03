#include "LibintCommunicator.hpp"

#include <libint2.hpp>


namespace libwint {


// Constructor
LibintCommunicator::LibintCommunicator() {
    libint2::initialize();
}

// Destructor: we don't want anyone else to possibly delete the singleton object
LibintCommunicator::~LibintCommunicator() {
    libint2::finalize();
}


/**
 *  @return the static singleton instance
 */
LibintCommunicator& LibintCommunicator::get() {
    static LibintCommunicator singleton_instance;  // Instantiated on first use and guaranteed to be destroyed
    return singleton_instance;
}


}  // namespace libwint
