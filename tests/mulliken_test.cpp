#define BOOST_TEST_MODULE "SOMullikenBasis_test"

#include "SOMullikenBasis.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( mulliken_test ) {
    Eigen::Matrix<double,7,7> aa;
    aa << 0.999999, 2.75445e-05, 2.10116e-16, -8.05074e-06, -2.79652e-18, 4.63697e-05, 3.64526e-17,
        2.75445e-05, 0.995582, -2.86904e-17, -0.00753316, 8.17833e-17,  0.00425522, -1.91957e-16,
        2.10116e-16, -2.86904e-17, 0.976792, 2.28438e-17 , 3.71927e-17,  3.65174e-16, -0.00047531,
        -8.05074e-06, -0.00753316, 2.28438e-17, 0.980476, 6.4507e-17, -0.0182836, 1.26353e-15,
        -2.79652e-18, 8.17833e-17, 3.71927e-17, 6.4507e-17, 0.999083, -4.46521e-15, -2.88329e-16,
        4.63697e-05, 0.00425522,  3.65174e-16, -0.0182836, -4.46521e-15,  0.0247316, -3.84334e-17,
        3.64526e-17, -1.91957e-16, -0.00047531, 1.26353e-15, -2.88329e-16, -3.84334e-17 , 0.0233358;
    Eigen::Matrix<double,7,7> C;
    C <<  -0.994435, -0.239159, 2.70425e-16, 0.0936834, 3.0127e-32, -0.11164, 1.19708e-15,
            -0.024097, 0.885736, -1.51665e-15, -0.479587, -4.1368e-31, 0.669578, -7.13724e-15,
            -7.16711e-19, -1.69007e-16, 0.607285, -2.60703e-15, 1.15279e-15, 9.63173e-15, 0.919233,
            -0.00316155, 0.0858963, 1.81567e-15, 0.74743, 2.31276e-31, 0.73849, -6.72651e-15,
            -2.04115e-34, 8.53707e-32, 1.68919e-15, -5.09334e-30, -1, -1.59442e-30, -1.14861e-16,
            0.00459374, 0.144039, 0.452998, 0.329472, 6.92404e-16, -0.709849, -0.732461,
            0.00459374, 0.144039, -0.452998, 0.329472, -6.92404e-16, -0.709849, 0.732461;

    Eigen::MatrixXd bb= aa;
    libwint::Molecule water ("../tests/ref_data/h2o.xyz");  // the relative path to the input .xyz-file w.r.t. the out-of-source build directory
    libwint::AOBasis ao_basis (water, "STO-3G");
    ao_basis.calculateIntegrals();
    libwint::SOMullikenBasis so_basis (ao_basis, C);
    std::cout<<std::endl;
    so_basis.calculateMullikenMatrix({0,1});
    std::cout<<so_basis.get_mulliken_matrix();

    // TO:DO add ref

}

