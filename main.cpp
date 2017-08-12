// First C++ project: interface LibInt2's calculations to matrices in eigen3

#include <libint2.hpp>


/**
 * Given an operator type, an orbital basis and atoms, calculates the one-body integrals (associated to that operator type)

 * @param opertype: a libint2::Operator (e.g. libint2::Operator::overlap)
 * @param obs:      a libint2::BasisSet object that represents the basis put on the molecule
 * @param atoms:    a std::vector<Atom>
 */
void compute_1body_integrals(const libint2::Operator &opertype, const libint2::BasisSet &obs, const std::vector<libint2::Atom> &atoms) {

    const auto nbf = obs.nbf();     // nbf: number of basis functions
    const auto nsh = obs.size();    // nsh: number of shells

    libint2::Engine engine (opertype, obs.max_nprim(), static_cast<int>(obs.max_l()));

    // Something extra for the nuclear attraction integrals
    if (opertype == libint2::Operator::nuclear) {
        engine.set_params(make_point_charges(atoms));
    }

    const auto shell2bf = obs.shell2bf();

    const auto& buffer = engine.results();      // vector that holds pointers to the shell sets
    // the values will change after every compute() call

    for(auto s1=0; s1!=nsh; ++s1) {
        for(auto s2=0; s2!=nsh; ++s2) {

            // std::cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
            engine.compute(obs[s1], obs[s2]);
            // std::cout << "done" << std::endl;
            auto ints_shellset = buffer[0];  // location of the computed integrals
            if (ints_shellset == nullptr)
                continue;  // nullptr returned if the entire shell-set was screened out

            auto bf1 = shell2bf[s1];  // first basis function in first shell
            auto n1 = obs[s1].size(); // number of basis functions in first shell
            auto bf2 = shell2bf[s2];  // first basis function in second shell
            auto n2 = obs[s2].size(); // number of basis functions in second shell

            // integrals are packed into ints_shellset in row-major (C) form
            // this iterates over integrals in this order
            for(auto f1=0; f1!=n1; ++f1) {
                for(auto f2=0; f2!=n2; ++f2) {
                    std::cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << std::endl;
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
    const auto xyzfilename = "/Users/laurentlemmens/Documents/Xcode/LibInt_CTF/LibInt_CTF/h2o.xyz";
    std::ifstream input_file(xyzfilename);
    auto atoms = libint2::read_dotxyz(input_file);
    libint2::BasisSet obs ("cc-pVDZ", atoms);




    // LEARNING AREA
    //print_shell_sizes(obs);

    compute_1body_integrals(libint2::Operator::overlap, obs, atoms);





    // Finalize libint2
    libint2::finalize();
    return 0;
}