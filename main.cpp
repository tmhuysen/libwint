#include <libint2.hpp>
#include <Eigen/Dense>      // <Eigen/Eigen> includes <Eigen/Dense> and <Eigen/Sparse>, so we might as well just include <Eigen/Dense> since we won't be using <Eigen/Sparse> in this code
#include <unsupported/Eigen/CXX11/Tensor>

#include "Molecule.hpp"
#include "Basis.hpp"


int main() {
    // Initialize libint2
    libint2::initialize();



    // 1. MOLECULE & BASIS SET SPECIFICATION

    const auto xyzfilename = "/Users/laurentlemmens/Software/libint-eigen/docs/h2o.xyz";
    std::string basis_name = "STO-3G";

    Molecule water (xyzfilename);


    // 2. CALCULATIONS
    Basis basis (water, basis_name);

    auto S = basis.compute_overlap_integrals();
    auto T = basis.compute_kinetic_integrals();
    auto V = basis.compute_nuclear_integrals();

    auto tei = basis.compute_two_electron_integrals();


    // Finalize libint2
    libint2::finalize();
    return 0;
}
