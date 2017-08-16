// First C++ project: interface LibInt2's calculations to matrices in eigen3

#include <libint2.hpp>


/**
 * Given an operator type, an orbital basis and atoms, calculates the one-body integrals (associated to that operator type)

 * @param opertype: a libint2::Operator (e.g. libint2::Operator::overlap)
 * @param obs:      a libint2::BasisSet object that represents the basis put on the molecule
 * @param atoms:    a std::vector<Atom>
 */
void compute_1body_integrals(const libint2::Operator &opertype, const libint2::BasisSet &obs, const std::vector<libint2::Atom> &atoms) {

    const auto nbf = obs.nbf();     // nbf: number of basis functions in the obs
    const auto nsh = obs.size();    // nsh: number of shells in the obs

    libint2::Engine engine (opertype, obs.max_nprim(), static_cast<int>(obs.max_l()));

    // Something extra for the nuclear attraction integrals
    if (opertype == libint2::Operator::nuclear) {
        engine.set_params(make_point_charges(atoms));
    }

    const auto shell2bf = obs.shell2bf();  // maps shell index to bf index

    const auto& buffer = engine.results();  // vector that holds pointers to computed shell sets
                                            // the values will change after every compute() call

    // One-body integrals are between two basis functions
    //  However, LibInt calculates integrals between shell sets, we will loop over the shells in the obs
    for(auto s1=0; s1!=nsh; ++s1) {
        for(auto s2=0; s2!=nsh; ++s2) {

            engine.compute(obs[s1], obs[s2]);
            auto ints_shellset = buffer[0];     // buffer is a vector that holds pointers to computed shell sets
                                                // if the zeroth element is nullptr, then the whole shell set has been exhausted
            if (ints_shellset == nullptr)
                continue;

            auto bf1 = shell2bf[s1];  // first basis function in first shell
            auto n1 = obs[s1].size(); // number of basis functions in first shell
            auto bf2 = shell2bf[s2];  // first basis function in second shell
            auto n2 = obs[s2].size(); // number of basis functions in second shell

            // integrals are packed into ints_shellset in row-major (C) form
            // this iterates over integrals in this order
            for(auto f1=0; f1!=n1; ++f1) {
                for(auto f2=0; f2!=n2; ++f2) {
                    auto computed_integral = ints_shellset[f1*n2 + f2];
                    std::cout << "  " << bf1+f1 << " " << bf2+f2 << " " << computed_integral << std::endl;
                }
            }
        }
    }
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


int main() {
    // Initialize libint2
    libint2::initialize();


    // 1. MOLECULE & BASIS SET SPECIFICATION
    const auto xyzfilename = "/Users/laurentlemmens/Software/LibInt_to_eigen3/docs/h2o.xyz";
    std::ifstream input_file(xyzfilename);
    auto atoms = libint2::read_dotxyz(input_file);
    libint2::BasisSet obs ("STO-3G", atoms);  // obs: orbital basis set




    // LEARNING AREA
    print_shell_sizes(obs);

    compute_1body_integrals(libint2::Operator::overlap, obs, atoms);





    // Finalize libint2
    libint2::finalize();
    return 0;
}
