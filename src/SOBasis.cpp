#include "SOBasis.hpp"

#include "transformations.hpp"



namespace libwint {

/*
 *  PRIVATE METHODS
 */

/**
 *  Parse a given FCIDUMP file for the one- and two-electron integrals
 */
void SOBasis::parseFCIDUMPFile(std::string fcidump_filename) {

    // Find the extension of the given path (https://stackoverflow.com/a/51992)
    std::string extension;
    std::string::size_type idx = fcidump_filename.rfind('.');

    if (idx != std::string::npos) {
        extension = fcidump_filename.substr(idx+1);
    } else {
        throw std::runtime_error("I did not find an extension in your given path.");
    }

    if (!(extension == "FCIDUMP")) {
        throw std::runtime_error("You did not provide a .FCIDUMP file name");
    }

    // If the xyz_filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
    std::ifstream input_file_stream (fcidump_filename);

    if (!input_file_stream.good()) {
        throw std::runtime_error("The provided FCIDUMP file is illegible. Maybe you specified a wrong path?");
    }



    // Do the actual parsing

    //  Get the number of orbitals to check if it's a valid FCIDUMP file
    std::string start_line;  // first line contains orbitals and electron count
    std::getline(input_file_stream, start_line);
    std::stringstream linestream (start_line);

    size_t value;
    char iter;

    while (linestream >> iter) {
        if (iter == '=') {
            linestream >> value;  // right here we have the number of orbitals
            if (this->K != value) {
                throw std::invalid_argument("It appears that the given number of spatial orbitals is inconsistent with the given FCIDUMP file.");
            }
            break;  // we can finish reading the linestream after we found K
        }
    }


    Eigen::MatrixXd h_SO = Eigen::MatrixXd::Zero(this->K, this->K);
    Eigen::Tensor<double, 4> g_SO (this->K, this->K, this->K, this->K);
    g_SO.setZero();

    //  Skip 3 lines
    for (size_t counter = 0; counter < 3; counter++) {
        std::getline(input_file_stream, start_line);
    }


    //  Start reading in the one- and two-electron integrals
    double x;
    size_t i, j, a, b;

    std::string line;
    while (std::getline(input_file_stream, line)) {
        std::istringstream iss (line);

        // Based on what the values of the indices are, we can read one-electron integrals, two-electron integrals and the internuclear repulsion energy
        //  See also (http://hande.readthedocs.io/en/latest/manual/integrals.html)
        iss >> x >> i >> a >> j >> b;

        //  Internuclear repulsion energy
        if ((i == 0) && (j == 0) && (a == 0) && (b == 0)) {
            this->internuclear_repulsion_energy = x;
        }

        //  Two-electron integrals are given in CHEMIST'S NOTATION
        if ((i > 0) && (a > 0) && (j > 0) && (b > 0)) {
            size_t p = i - 1;
            size_t q = a - 1;
            size_t r = j - 1;
            size_t s = b - 1;
            g_SO(p,q,r,s) = x;

            // Set the integral indices for the permutational symmetries as well
            g_SO(p,q,s,r) = x;
            g_SO(q,p,r,s) = x;
            g_SO(q,p,s,r) = x;

            g_SO(r,s,p,q) = x;
            g_SO(s,r,p,q) = x;
            g_SO(r,s,q,p) = x;
            g_SO(s,r,q,p) = x;
        }

        //  One-electron integrals
        if ((i > 0) && (a > 0) && (j == 0)) {
            size_t p = i - 1;
            size_t q = a - 1;
            h_SO(p,q) = x;

        }


    }

}


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param atomic orbital instance and a coefficient matrix @param C (i.e. a basis transformation matrix) that links the SO basis to the AO basis
 */
SOBasis::SOBasis(const libwint::AOBasis& ao_basis, const Eigen::MatrixXd& C) :
        K(ao_basis.calculateNumberOfBasisFunctions())
{

    this->internuclear_repulsion_energy = ao_basis.get_molecule().calculateInternuclearRepulsionEnergy();


    Eigen::MatrixXd h_AO = ao_basis.get_T() + ao_basis.get_V();
    this->h_SO = libwint::transformations::transform_AO_to_SO(h_AO, C);

    this->g_SO = libwint::transformations::transform_AO_to_SO(ao_basis.get_g(), C);
}

/**
 *  Constructor based on a given path to an FCIDUMP file
 */
SOBasis::SOBasis(std::string fcidump_filename, size_t K) :
        K(K)
{

}


/*
 *  PUBLIC METHODS
 */

/**
 *  Transform the one- and two-electron integrals according to the basis transformation matrix @param T
 */
void SOBasis::transform(const Eigen::MatrixXd& T) {

    this->h_SO = libwint::transformations::transformOneElectronIntegrals(this->h_SO, T);
    this->g_SO = libwint::transformations::transformTwoElectronIntegrals(this->g_SO, T);
}


/**
 *  Transform the one- and two-electron integrals according to the Jacobi rotation parameters p, q and a given angle theta in radians.
 */
void SOBasis::rotateJacobi(size_t p, size_t q, double theta) {

    // We can use our specialized rotate{One,Two}ElectronIntegralsJacobi functions
    this->h_SO = libwint::transformations::rotateOneElectronIntegralsJacobi(this->h_SO, p, q, theta);
    this->g_SO = libwint::transformations::rotateTwoElectronIntegralsJacobi(this->g_SO, p, q, theta);
}


}  // namespace libwint
