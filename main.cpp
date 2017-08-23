#include <libint2.hpp>
#include <Eigen/Dense>      // <Eigen/Eigen> includes <Eigen/Dense> and <Eigen/Sparse>, so we might as well just include <Eigen/Dense> since we won't be using <Eigen/Sparse> in this code
#include <unsupported/Eigen/CXX11/Tensor>

#include "libint-wrapper.hpp"
#include "Molecule.hpp"

int main() {
    // Initialize libint2
    libint2::initialize();


    // 1. MOLECULE & BASIS SET SPECIFICATION

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";

    Molecule water (xyzfilename);
    for(auto atom : water.atoms) {
        std::cout << atom.atomic_number << std::endl;
    }



    // std::ifstream input_file (xyzfilename);
    // auto atoms = libint2::read_dotxyz(input_file);
    // libint2::BasisSet obs ("STO-3G", atoms);    // obs: orbital basis set
                                                // a libint2::BasisSet is a decorated std::vector<libint2::Shell>

    // std::cout << obs.size() << std::endl;
    // std::cout << obs.nbf() << std::endl;

    // 2. CALCULATE ONE- AND TWO BODY INTEGRALS
    //auto S = compute_1body_integrals(libint2::Operator::overlap, obs, atoms);
    // auto T = compute_1body_integrals(libint2::Operator::kinetic, obs, atoms);
    // auto V = compute_1body_integrals(libint2::Operator::nuclear, obs, atoms);

    // std::cout << S << std::endl << T << std::endl << V << std::endl;

    // auto tei = compute_2body_integrals(obs, atoms);

    // compute_2body_integrals(obs, atoms);

    // Finalize libint2
    libint2::finalize();
    return 0;
}
