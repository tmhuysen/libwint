//
// Created by Laurent Lemmens on 21/08/17.
//

#include <libint2.hpp>
#include <Eigen/Dense>      // <Eigen/Eigen> includes <Eigen/Dense> and <Eigen/Sparse>, so we might as well just include <Eigen/Dense> since we won't be using <Eigen/Sparse> in this code
#include "libint-wrapper.hpp"

/**
 * Given an operator type, an orbital basis and atoms, calculates the one-body integrals (associated to that operator type)

 * @param opertype: a libint2::Operator (e.g. libint2::Operator::overlap)
 * @param obs:      a libint2::BasisSet object that represents the basis put on the molecule
 * @param atoms:    a std::vector<Atom>
 */
Eigen::MatrixXf compute_1body_integrals(const libint2::Operator &opertype, const libint2::BasisSet &obs, const std::vector<libint2::Atom> &atoms) {

    const auto nsh = obs.size();    // nsh: number of shells in the obs
    const auto nbf = obs.nbf();     // nbf: number of basis functions in the obs

    // Initialize the eigen matrix:
    //  Since the matrices we will encounter (S, T, V) are symmetric, the issue of row major vs column major doesn't matter.
    Eigen::MatrixXf M_result (nbf, nbf);

    // Construct the libint2 engine
    libint2::Engine engine(opertype, obs.max_nprim(), static_cast<int>(obs.max_l()));
    //  Something extra for the nuclear attraction integrals
    if (opertype == libint2::Operator::nuclear) {
        engine.set_params(make_point_charges(atoms));
    }

    const auto shell2bf = obs.shell2bf();  // maps shell index to bf index

    const auto &buffer = engine.results();  // vector that holds pointers to computed shell sets
    // actually, buffer.size() is always 1, so buffer[0] is a pointer to
    //      the first calculated integral of these specific shells
    // the values that buffer[0] points to will change after every compute() call

    // One-body integrals are between two basis functions, so we'll need two loops.
    // However, LibInt calculates integrals between libint2::Shells, we will loop over the shells (sh) in the obs
    for (auto sh1 = 0; sh1 != nsh; ++sh1) {  // sh1: shell 1
        for (auto sh2 = 0; sh2 != nsh; ++sh2) {  // sh2: shell2
            // Calculate integrals between the two shells (obs is a decorated std::vector<libint2::Shell>)
            engine.compute(obs[sh1], obs[sh2]);

            auto calculated_integrals = buffer[0];  // is actually a pointer: const double *

            if (calculated_integrals == nullptr)    // if the zeroth element is nullptr, then the whole shell has been exhausted
                continue;


            // Extract the calculated integrals from calculated_integrals.
            // In calculated_integrals, the integrals are stored in row major form.
            auto bf1 = shell2bf[sh1];  // (index of) first bf in sh1
            auto bf2 = shell2bf[sh2];  // (index of) first bf in sh2

            auto nbf_sh1 = obs[sh1].size();  // number of basis functions in first shell
            auto nbf_sh2 = obs[sh2].size();  // number of basis functions in second shell

            for (auto f1 = 0; f1 != nbf_sh1; ++f1) {     // f1: index of basis function within shell 1
                for (auto f2 = 0; f2 != nbf_sh2; ++f2) { // f2: index of basis function within shell 2
                    auto computed_integral = calculated_integrals[f1 * nbf_sh2 + f2];  // integrals are packed in row-major form
                    M_result(bf1 + f1, bf2 + f2) = computed_integral;
                }
            }

        }
    }

    return M_result;
}



/**
 * Prints the sizes (i.e. the number of basis functions in them) of all shells in a given basis set object.
 *
 * @param obs:  the given basis set object
 */
void print_shell_sizes(const libint2::BasisSet &obs) {
    auto size = obs.size();
    for (auto i = 0; i < size; i++) {
        auto sh = obs[i];               // i traverses the shell
        std::cout << "Shell nr.: " << i << "\t Shell size: " << sh.size() << std::endl;
    }
}


void test(int a){
    std::cout << a << std::endl;
}
