#include <fstream>
#include <stdlib.h>  /* getenv */

#include "gtest/gtest.h"

#include "xcint.h"
#include "numgrid.h"

TEST(xcint, energy_spherical)
{
    int ierr;
    double dot;

    int num_centers = 2;

    double *center_coordinates = NULL;
    center_coordinates = new double[3*num_centers];
    center_coordinates[0] = 1.7;
    center_coordinates[1] = 0.0;
    center_coordinates[2] = 0.0;
    center_coordinates[3] = 0.0;
    center_coordinates[4] = 0.0;
    center_coordinates[5] = 0.0;

    int *center_elements = NULL;
    center_elements = new int[num_centers];
    center_elements[0] = 9;
    center_elements[1] = 1;

    int num_shells = 9;

    int *shell_centers = NULL;
    shell_centers = new int[num_shells];
    int *shell_l_quantum_numbers = NULL;
    shell_l_quantum_numbers = new int[num_shells];
    int *shell_num_primitives = NULL;
    shell_num_primitives = new int[num_shells];
    shell_centers[0] = 1;
    shell_l_quantum_numbers[0] = 0;
    shell_num_primitives[0] = 9;
    shell_centers[1] = 1;
    shell_l_quantum_numbers[1] = 0;
    shell_num_primitives[1] = 9;
    shell_centers[2] = 1;
    shell_l_quantum_numbers[2] = 0;
    shell_num_primitives[2] = 1;
    shell_centers[3] = 1;
    shell_l_quantum_numbers[3] = 1;
    shell_num_primitives[3] = 4;
    shell_centers[4] = 1;
    shell_l_quantum_numbers[4] = 1;
    shell_num_primitives[4] = 1;
    shell_centers[5] = 1;
    shell_l_quantum_numbers[5] = 2;
    shell_num_primitives[5] = 1;
    shell_centers[6] = 2;
    shell_l_quantum_numbers[6] = 0;
    shell_num_primitives[6] = 4;
    shell_centers[7] = 2;
    shell_l_quantum_numbers[7] = 0;
    shell_num_primitives[7] = 1;
    shell_centers[8] = 2;
    shell_l_quantum_numbers[8] = 1;
    shell_num_primitives[8] = 1;

    int num_exponents = 0;
    for (int i = 0; i < num_shells; i++)
    {
        num_exponents += shell_num_primitives[i];
    }

    double *primitive_exponents = NULL;
    primitive_exponents = new double[num_exponents];
    double *contraction_coefficients = NULL;
    contraction_coefficients = new double[num_exponents];
    primitive_exponents[0] =   1.471000000000e+04;
    contraction_coefficients[0] =   6.863650000000e-01;
    primitive_exponents[1] =   2.207000000000e+03;
    contraction_coefficients[1] =   1.274350000000e+00;
    primitive_exponents[2] =   5.028000000000e+02;
    contraction_coefficients[2] =   2.139130000000e+00;
    primitive_exponents[3] =   1.426000000000e+02;
    contraction_coefficients[3] =   3.130550000000e+00;
    primitive_exponents[4] =   4.647000000000e+01;
    contraction_coefficients[4] =   3.638230000000e+00;
    primitive_exponents[5] =   1.670000000000e+01;
    contraction_coefficients[5] =   2.641480000000e+00;
    primitive_exponents[6] =   6.356000000000e+00;
    contraction_coefficients[6] =   7.553570000000e-01;
    primitive_exponents[7] =   1.316000000000e+00;
    contraction_coefficients[7] =   1.342700000000e-02;
    primitive_exponents[8] =   3.897000000000e-01;
    contraction_coefficients[8] =  -8.197600000000e-04;
    primitive_exponents[9] =   1.471000000000e+04;
    contraction_coefficients[9] =  -1.570740000000e-01;
    primitive_exponents[10] =   2.207000000000e+03;
    contraction_coefficients[10] =  -3.001720000000e-01;
    primitive_exponents[11] =   5.028000000000e+02;
    contraction_coefficients[11] =  -4.915140000000e-01;
    primitive_exponents[12] =   1.426000000000e+02;
    contraction_coefficients[12] =  -7.849910000000e-01;
    primitive_exponents[13] =   4.647000000000e+01;
    contraction_coefficients[13] =  -9.347560000000e-01;
    primitive_exponents[14] =   1.670000000000e+01;
    contraction_coefficients[14] =  -1.005480000000e+00;
    primitive_exponents[15] =   6.356000000000e+00;
    contraction_coefficients[15] =  -3.204660000000e-01;
    primitive_exponents[16] =   1.316000000000e+00;
    contraction_coefficients[16] =   4.928530000000e-01;
    primitive_exponents[17] =   3.897000000000e-01;
    contraction_coefficients[17] =   1.999410000000e-01;
    primitive_exponents[18] =   3.897000000000e-01;
    contraction_coefficients[18] =   3.515260000000e-01;
    primitive_exponents[19] =   2.267000000000e+01;
    contraction_coefficients[19] =   3.164380000000e+00;
    primitive_exponents[20] =   4.977000000000e+00;
    contraction_coefficients[20] =   2.497710000000e+00;
    primitive_exponents[21] =   1.347000000000e+00;
    contraction_coefficients[21] =   1.051860000000e+00;
    primitive_exponents[22] =   3.471000000000e-01;
    contraction_coefficients[22] =   1.739750000000e-01;
    primitive_exponents[23] =   3.471000000000e-01;
    contraction_coefficients[23] =   3.797590000000e-01;
    primitive_exponents[24] =   1.640000000000e+00;
    contraction_coefficients[24] =   6.775590000000e+00;
    primitive_exponents[25] =   1.301000000000e+01;
    contraction_coefficients[25] =   9.610660000000e-02;
    primitive_exponents[26] =   1.962000000000e+00;
    contraction_coefficients[26] =   1.630200000000e-01;
    primitive_exponents[27] =   4.446000000000e-01;
    contraction_coefficients[27] =   1.855450000000e-01;
    primitive_exponents[28] =   1.220000000000e-01;
    contraction_coefficients[28] =   7.374380000000e-02;
    primitive_exponents[29] =   1.220000000000e-01;
    contraction_coefficients[29] =   1.471230000000e-01;
    primitive_exponents[30] =   7.270000000000e-01;
    contraction_coefficients[30] =   9.568810000000e-01;

    // generate grid
    context_t *numgrid_context = new_context();
    double radial_precision = 1.0e-12;
    int min_num_angular_points = 86;
    int max_num_angular_points = 302;
    int num_outer_centers = 0;
    double *outer_center_coordinates = NULL;
    int *outer_center_elements = NULL;
    ierr = generate_grid(numgrid_context,
                         radial_precision,
                         min_num_angular_points,
                         max_num_angular_points,
                         num_centers,
                         center_coordinates,
                         center_elements,
                         num_outer_centers,
                         outer_center_coordinates,
                         outer_center_elements,
                         num_shells,
                         shell_centers,
                         shell_l_quantum_numbers,
                         shell_num_primitives,
                         primitive_exponents);
    int num_points = get_num_points(numgrid_context);
    ASSERT_EQ(num_points, 31424);

    double *grid = (double*) get_grid(numgrid_context);

    xcint_context_t *xcint_context = xcint_new();

    // we do this twice to test idempotent xcint_set_basis
    ierr = xcint_set_basis(xcint_context,
                           XCINT_BASIS_SPHERICAL,
                           num_centers,
                           center_coordinates,
                           num_shells - 1, // wrong on purpose
                           shell_centers,
                           shell_l_quantum_numbers,
                           shell_num_primitives,
                           primitive_exponents,
                           contraction_coefficients);

    ierr = xcint_set_basis(xcint_context,
                           XCINT_BASIS_SPHERICAL,
                           num_centers,
                           center_coordinates,
                           num_shells,
                           shell_centers,
                           shell_l_quantum_numbers,
                           shell_num_primitives,
                           primitive_exponents,
                           contraction_coefficients);

    delete[] center_coordinates;
    center_coordinates = NULL;
    delete[] center_elements;
    center_elements = NULL;
    delete[] shell_centers;
    shell_centers = NULL;
    delete[] shell_l_quantum_numbers;
    shell_l_quantum_numbers = NULL;
    delete[] shell_num_primitives;
    shell_num_primitives = NULL;
    delete[] primitive_exponents;
    primitive_exponents = NULL;
    delete[] contraction_coefficients;
    contraction_coefficients = NULL;

    int mat_dim = 19;

    double *dmat = NULL;
    dmat = new double[mat_dim*mat_dim];
    std::fill(&dmat[0], &dmat[mat_dim*mat_dim], 0.0);

    char* test_directory = getenv("XCINT_TEST_DIRECTORY");
    if (test_directory == NULL)
    {
        fputs("ERROR: env variable XCINT_TEST_DIRECTORY not set!\n", stderr);
        abort();
    }
    std::string dmat_file_name = test_directory + std::string("/dmat.txt");
    std::ifstream infile(dmat_file_name.c_str());
    int i;
    double d;
    while (infile >> i >> d)
    {
        dmat[i] = d;
    }

    double *vxc = NULL;
    vxc = new double[mat_dim*mat_dim];

    // we call it twice to test idempotency
    ierr = xcint_set_functional(xcint_context, "lda");
    ierr = xcint_set_functional(xcint_context, "lda");

    double exc = 0.0;
    double num_electrons = 0.0;

    ierr = xcint_integrate_scf(xcint_context,
                               XCINT_MODE_RKS,
                               num_points,
                               grid,
                               dmat,
                               &exc,
                               vxc,
                               &num_electrons);

    ASSERT_NEAR(num_electrons, 9.999992074832, 1.0e-11);
    ASSERT_NEAR(exc, -20.421064966255539, 1.0e-11);

    dot = 0.0;
    for (int i = 0; i < mat_dim*mat_dim; i++)
    {
        dot += vxc[i]*dmat[i];
    }
    ASSERT_NEAR(dot, -6.729996811122, 1.0e-11);

    ierr = xcint_set_functional(xcint_context, "b3lyp");

    ierr = xcint_integrate_scf(xcint_context,
                               XCINT_MODE_RKS,
                               num_points,
                               grid,
                               dmat,
                               &exc,
                               vxc,
                               &num_electrons);

    ASSERT_NEAR(num_electrons, 9.999992074832, 1.0e-11);
    ASSERT_NEAR(exc, -17.475254754458547, 1.0e-11);

    dot = 0.0;
    for (int i = 0; i < mat_dim*mat_dim; i++)
    {
        dot += vxc[i]*dmat[i];
    }
    ASSERT_NEAR(dot, -5.6105711653099748, 1.0e-11);

    delete[] dmat;
    dmat = NULL;
    delete[] vxc;
    vxc = NULL;

    free_context(numgrid_context);
    xcint_free(xcint_context);
}
