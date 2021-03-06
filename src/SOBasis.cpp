#include "SOBasis.hpp"




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
        //  I think the documentation is a bit unclear for the two-electron integrals, but we can rest assured that FCIDUMP files give the two-electron integrals in CHEMIST's notation.
        iss >> x >> i >> a >> j >> b;

        //  Internuclear repulsion energy (skipped)
        if ((i == 0) && (j == 0) && (a == 0) && (b == 0)) {}

        //  Single-particle eigenvalues (skipped)
        else if ((a == 0) && (j == 0) && (b == 0)) {}

        //  One-electron integrals (h_core)
        else if ((j == 0) && (b == 0)) {
            size_t p = i - 1;
            size_t q = a - 1;
            h_SO(p,q) = x;

            // Apply the permutational symmetry for real orbitals
            h_SO(q,p) = x;
        }

        //  Two-electron integrals are given in CHEMIST'S NOTATION, so just copy them over
        else if ((i > 0) && (a > 0) && (j > 0) && (b > 0)) {
            size_t p = i - 1;
            size_t q = a - 1;
            size_t r = j - 1;
            size_t s = b - 1;
            g_SO(p,q,r,s) = x;

            // Apply the permutational symmetries for real orbitals
            g_SO(p,q,s,r) = x;
            g_SO(q,p,r,s) = x;
            g_SO(q,p,s,r) = x;

            g_SO(r,s,p,q) = x;
            g_SO(s,r,p,q) = x;
            g_SO(r,s,q,p) = x;
            g_SO(s,r,q,p) = x;
        }
    }  // while loop

    this->h_SO = h_SO;
    this->g_SO = g_SO;
}


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param atomic orbital instance and a coefficient matrix @param C (i.e. a basis transformation matrix) that links the SO basis to the AO basis
 */
SOBasis::SOBasis(const libwint::AOBasis& ao_basis, const Eigen::MatrixXd& C) :
        K (ao_basis.calculateNumberOfBasisFunctions())
{

    Eigen::MatrixXd h_AO = ao_basis.get_T() + ao_basis.get_V();
    this->h_SO = libwint::transformations::transform_AO_to_SO(h_AO, C);

    this->g_SO = libwint::transformations::transform_AO_to_SO(ao_basis.get_g(), C);
}

/**
 *  Constructor based on a given path to an FCIDUMP file
 */
SOBasis::SOBasis(std::string fcidump_filename, size_t K, bool hack) :
        K (K)
{

    // This sets this->h_SO, this->g_SO, this->internuclear_repulsion by reading in the FCIDUMP file
    if(hack){
        this->parseFCIDUMPFile(fcidump_filename);
    }
    else{
        this->parseOne(fcidump_filename);
        this->parseTwo(fcidump_filename);
    }
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

/**
 *  Parse a given One file for the one- and two-electron integrals and overlap.
 */
void SOBasis::parseTwo(std::string fcidump_filename) {
    std::ifstream input_file_stream (fcidump_filename + ".two");
    if (!input_file_stream.good()) {
        throw std::runtime_error("The provided BLANKKKKKK file is illegible. Maybe you specified a wrong path?");
    }
    //  Start reading in the one- and two-electron integrals
    std::string dtoe;
    size_t i, j, a, b;
    double x;
    Eigen::Tensor<double, 4> g_SO (this->K, this->K, this->K, this->K);
    g_SO.setZero();
    std::string line;
    while (std::getline(input_file_stream, line)) {
        std::istringstream iss(line);

        // Based on what the values of the indices are, we can read one-electron integrals, two-electron integrals and the internuclear repulsion energy
        //  See also (http://hande.readthedocs.io/en/latest/manual/integrals.html)
        //  I think the documentation is a bit unclear for the two-electron integrals, but we can rest assured that FCIDUMP files give the two-electron integrals in CHEMIST's notation.
        iss >> i >> a >> j >> b >> dtoe;
        //std::cout << dtoe;
        auto e(dtoe.find_first_of("Dd"));
        if (e != std::string::npos) {
            dtoe[e] = 'E';
        }
        x = std::stod(dtoe);
        //  Two-electron integrals are given in CHEMIST'S NOTATION, so just copy them over
        if ((i > 0) && (a > 0) && (j > 0) && (b > 0)) {
            size_t p = i - 1;
            size_t q = a - 1;
            size_t r = j - 1;
            size_t s = b - 1;
            g_SO(p,q,r,s) = x;

        }
    }  // while loop
    this->g_SO = g_SO;
}

void SOBasis::parseOne(std::string fcidump_filename) {


    std::ifstream input_file_stream (fcidump_filename + ".one");
    if (!input_file_stream.good()) {
        throw std::runtime_error("The provided BLANKKKKKK file is illegible. Maybe you specified a wrong path?");
    }
    //  Start reading in the one- and two-electron integrals
    std::string dtoe;
    size_t i, j;
    double x;
    Eigen::MatrixXd h_SO = Eigen::MatrixXd::Zero(this->K, this->K);
    std::string line;
    while (std::getline(input_file_stream, line)) {
        std::istringstream iss(line);

        // Based on what the values of the indices are, we can read one-electron integrals, two-electron integrals and the internuclear repulsion energy
        //  See also (http://hande.readthedocs.io/en/latest/manual/integrals.html)
        //  I think the documentation is a bit unclear for the two-electron integrals, but we can rest assured that FCIDUMP files give the two-electron integrals in CHEMIST's notation.
        iss >> i >> j >> dtoe;
        auto e(dtoe.find_first_of("Dd"));
        if (e != std::string::npos) {
            dtoe[e] = 'E';
        }
        x = std::stod(dtoe);
        //  Two-electron integrals are given in CHEMIST'S NOTATION, so just copy them over
        if ((i > 0) && (j > 0)) {
            size_t p = i - 1;
            size_t r = j - 1;
            h_SO(p,r) = x;

        }
    }  // while loop
    this->h_SO = h_SO;
}


}  // namespace libwint
