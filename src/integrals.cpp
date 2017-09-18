//
// Created by Laurent Lemmens on 21/08/17.
//

#include "integrals.hpp"

/**
 * Given an operator type, an orbital basis and atoms, calculates the one-body integrals (associated to that operator type)

 * @param opertype: a libint2::Operator (e.g. libint2::Operator::overlap)
 * @param obs:      a libint2::BasisSet object that represents the basis put on the molecule
 * @param atoms:    a std::vector<Atom>

 * @return: an Eigen::MatrixXd storing the integrals
 */
Eigen::MatrixXd compute_1body_integrals(const libint2::Operator& opertype, const libint2::BasisSet& obs, const std::vector<libint2::Atom>& atoms) {

    const auto nsh = obs.size();    // nsh: number of shells in the obs
    const auto nbf = obs.nbf();     // nbf: number of basis functions in the obs

    // Initialize the eigen matrix:
    //  Since the matrices we will encounter (S, T, V) are symmetric, the issue of row major vs column major doesn't matter.
    Eigen::MatrixXd M_result (nbf, nbf);

    // Construct the libint2 engine
    libint2::Engine engine (opertype, obs.max_nprim(), static_cast<int>(obs.max_l()));
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
        for (auto sh2 = 0; sh2 != nsh; ++sh2) {  // sh2: shell 2
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
                    auto computed_integral = calculated_integrals[f2 + f1 * nbf_sh2];  // integrals are packed in row-major form
                    M_result(bf1 + f1, bf2 + f2) = computed_integral;
                }
            }

        }
    }

    return M_result;
}


/**
 * Calculates the two-electron integrals, given an orbital basis and atoms.

 * @param obs:      a libint2::BasisSet object that represents the basis put on the molecule
 * @param atoms:    a std::vector<Atom>

 * @return: an Eigen::Tensor<double, 4> storing the integrals
 */
Eigen::Tensor<double, 4> compute_2body_integrals(const libint2::BasisSet& obs, const std::vector<libint2::Atom>& atoms){
    // We have to static_cast to LONG, as clang++ else gives the following errors:
    //  error: non-constant-expression cannot be narrowed from type 'unsigned long' to 'value_type' (aka 'long') in initializer list
    //  note: insert an explicit cast to silence this issue

    const auto nsh = static_cast<long>(obs.size());
    const auto nbf = obs.nbf();

    // Initialize the two-electron integrals tensor
    Eigen::Tensor<double, 4> tei (nbf, nbf, nbf, nbf);       // Create a rank 4 tensor with dimensions nbf


    // Construct the libint2 engine
    libint2::Engine engine (libint2::Operator::coulomb, obs.max_nprim(), static_cast<int>(obs.max_l()));

    const auto shell2bf = obs.shell2bf();  // maps shell index to bf index

    const auto &buffer = engine.results();  // vector that holds pointers to computed shell sets
    // actually, buffer.size() is always 1, so buffer[0] is a pointer to
    //      the first calculated integral of these specific shells
    // the values that buffer[0] points to will change after every compute() call

    // Two-electron integrals are between four basis functions, so we'll need four loops.
    // However, LibInt calculates integrals between libint2::Shells, we will loop over the shells (sh) in the obs
    for (auto sh1 = 0; sh1 != nsh; ++sh1) {  // sh1: shell 1
        for (auto sh2 = 0; sh2 != nsh; ++sh2) {  // sh2: shell 2
            for (auto sh3=0; sh3 != nsh; ++sh3) {  // sh3: shell 3
                for (auto sh4=0; sh4 != nsh; ++sh4) {  //sh4: shell 4
                    // Calculate integrals between the two shells (obs is a decorated std::vector<libint2::Shell>)
                    engine.compute(obs[sh1], obs[sh2], obs[sh3], obs[sh4]);

                    auto calculated_integrals = buffer[0];

                    if (calculated_integrals == nullptr)    // if the zeroth element is nullptr, then the whole shell has been exhausted
                        continue;

                    // Extract the calculated integrals from calculated_integrals.
                    // In calculated_integrals, the integrals are stored in row major form.
                    auto bf1 = static_cast<long>(shell2bf[sh1]);  // (index of) first bf in sh1
                    auto bf2 = static_cast<long>(shell2bf[sh2]);  // (index of) first bf in sh2
                    auto bf3 = static_cast<long>(shell2bf[sh3]);  // (index of) first bf in sh3
                    auto bf4 = static_cast<long>(shell2bf[sh4]);  // (index of) first bf in sh4


                    auto nbf_sh1 = static_cast<long>(obs[sh1].size());  // number of basis functions in first shell
                    auto nbf_sh2 = static_cast<long>(obs[sh2].size());  // number of basis functions in second shell
                    auto nbf_sh3 = static_cast<long>(obs[sh3].size());  // number of basis functions in third shell
                    auto nbf_sh4 = static_cast<long>(obs[sh4].size());  // number of basis functions in fourth shell

                    for (auto f1 = 0L; f1 != nbf_sh1; ++f1) {
                        for (auto f2 = 0L; f2 != nbf_sh2; ++f2) {
                            for (auto f3 = 0L; f3 != nbf_sh3; ++f3) {
                                for (auto f4 = 0L; f4 != nbf_sh4; ++f4) {
                                    auto computed_integral = calculated_integrals[f4 + nbf_sh4 * (f3 + nbf_sh3 * (f2 + nbf_sh2 * (f1)))];  // row-major storage accessing
                                    tei(f1+bf1, f2+bf2, f3+bf3, f4+bf4) = computed_integral;
                                }
                            }
                        }
                    }


                }
            }
        }
    }

    return tei;
};

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
