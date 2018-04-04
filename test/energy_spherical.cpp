#include <fstream>
#include <cstdlib>  /* getenv */

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

    double x_coordinates_bohr[2];
    x_coordinates_bohr[0] = 1.7;
    x_coordinates_bohr[1] = 0.0;

    double y_coordinates_bohr[2];
    y_coordinates_bohr[0] = 0.0;
    y_coordinates_bohr[1] = 0.0;

    double z_coordinates_bohr[2];
    z_coordinates_bohr[0] = 0.0;
    z_coordinates_bohr[1] = 0.0;

    int *proton_charges = NULL;
    proton_charges = new int[num_centers];
    proton_charges[0] = 9;
    proton_charges[1] = 1;

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

    // FH molecule in cc-pVDZ basis
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

    double alpha_max[2];
    alpha_max[0] = 14710.0;
    alpha_max[1] = 13.01;

    double alpha_min[2][3];
    alpha_min[0][0] = 0.3897;
    alpha_min[0][1] = 0.3471;
    alpha_min[0][2] = 1.64;
    alpha_min[1][0] = 0.122;
    alpha_min[1][1] = 0.727;
    alpha_min[1][2] = 0.0; // not used

    // generate grid
    double radial_precision = 1.0e-12;
    int min_num_angular_points = 86;
    int max_num_angular_points = 302;

    int max_l_quantum_numbers[2];
    max_l_quantum_numbers[0] = 2;
    max_l_quantum_numbers[1] = 1;

    int num_points = 31424;  // cheat to avoid reading twice
    double *grid_x_bohr = new double[num_points];
    double *grid_y_bohr = new double[num_points];
    double *grid_z_bohr = new double[num_points];
    double *grid_w = new double[num_points];
    num_points = 0;

    for (int center_index = 0; center_index < num_centers; center_index++)
    {
        context_t *context =
            numgrid_new_atom_grid(radial_precision,
                                  min_num_angular_points,
                                  max_num_angular_points,
                                  proton_charges[center_index],
                                  alpha_max[center_index],
                                  max_l_quantum_numbers[center_index],
                                  alpha_min[center_index]);

        int num_points_center = numgrid_get_num_grid_points(context);

        double *atom_grid_x_bohr = new double[num_points_center];
        double *atom_grid_y_bohr = new double[num_points_center];
        double *atom_grid_z_bohr = new double[num_points_center];
        double *atom_grid_w = new double[num_points_center];

        numgrid_get_grid(context,
                         num_centers,
                         center_index,
                         x_coordinates_bohr,
                         y_coordinates_bohr,
                         z_coordinates_bohr,
                         proton_charges,
                         atom_grid_x_bohr,
                         atom_grid_y_bohr,
                         atom_grid_z_bohr,
                         atom_grid_w);

        for (int i = 0; i < num_points_center; i++)
        {
            grid_x_bohr[num_points + i] = atom_grid_x_bohr[i];
            grid_y_bohr[num_points + i] = atom_grid_y_bohr[i];
            grid_z_bohr[num_points + i] = atom_grid_z_bohr[i];
            grid_w[num_points + i] = atom_grid_w[i];
        }
        num_points += num_points_center;

        delete[] atom_grid_x_bohr;
        delete[] atom_grid_y_bohr;
        delete[] atom_grid_z_bohr;
        delete[] atom_grid_w;

        numgrid_free_atom_grid(context);
    }

    ASSERT_EQ(num_points, 31424);

    xcint_context_t *xcint_context = xcint_new_context();

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
    delete[] proton_charges;
    proton_charges = NULL;
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
                               grid_x_bohr,
                               grid_y_bohr,
                               grid_z_bohr,
                               grid_w,
                               dmat,
                               &exc,
                               vxc,
                               &num_electrons);

    ASSERT_NEAR(num_electrons, 9.999992072209077, 1.0e-12);
    ASSERT_NEAR(exc, -20.421064966253642, 1.0e-12);

    dot = 0.0;
    for (int i = 0; i < mat_dim*mat_dim; i++)
    {
        dot += vxc[i]*dmat[i];
    }
    ASSERT_NEAR(dot, -6.729996811121003, 1.0e-12);

    ierr = xcint_set_functional(xcint_context, "b3lyp");

    ierr = xcint_integrate_scf(xcint_context,
                               XCINT_MODE_RKS,
                               num_points,
                               grid_x_bohr,
                               grid_y_bohr,
                               grid_z_bohr,
                               grid_w,
                               dmat,
                               &exc,
                               vxc,
                               &num_electrons);

    ASSERT_NEAR(num_electrons, 9.999992072209077, 1.0e-12);
    ASSERT_NEAR(exc, -17.475254754225027, 1.0e-12);

    dot = 0.0;
    for (int i = 0; i < mat_dim*mat_dim; i++)
    {
        dot += vxc[i]*dmat[i];
    }
    ASSERT_NEAR(dot, -5.610571165249672, 1.0e-12);

    delete[] dmat;
    dmat = NULL;
    delete[] vxc;
    vxc = NULL;

    xcint_free_context(xcint_context);
}
