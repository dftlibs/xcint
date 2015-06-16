#include <cmath>

#include "XCint.h"
#include "xcint.h"
#include "xcint_c_parameters.h"
#include "MemAllocator.h"
#include "numgrid.h"

int main(int argc, char** argv)
{
    int return_code = 0;
    int ierr;

    xcint_context_t *xcint_context = xcint_new();

    int num_centers;
    int num_shells;
    double *center_coordinates = NULL;
    int *center_elements = NULL;
    int *shell_centers = NULL;
    int *shell_l_quantum_numbers = NULL;
    int *shell_num_primitives = NULL;
    double *primitive_exponents = NULL;
    double *contraction_coefficients = NULL;
    size_t block_size;
    num_centers = 2;
    num_shells = 9;
    block_size = 3*num_centers*sizeof(double);
    center_coordinates = (double*) MemAllocator::allocate(block_size);
    center_coordinates[0] =   1.700000000000e+00;
    center_coordinates[1] =   0.000000000000e+00;
    center_coordinates[2] =   0.000000000000e+00;
    center_coordinates[3] =   0.000000000000e+00;
    center_coordinates[4] =   0.000000000000e+00;
    center_coordinates[5] =   0.000000000000e+00;
    block_size = num_centers*sizeof(int);
    center_elements = (int*) MemAllocator::allocate(block_size);
    center_elements[0] = 9;
    center_elements[1] = 1;
    block_size = num_shells*sizeof(int);
    shell_centers = (int*) MemAllocator::allocate(block_size);
    shell_l_quantum_numbers = (int*) MemAllocator::allocate(block_size);
    shell_num_primitives = (int*) MemAllocator::allocate(block_size);
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
    block_size = 31*sizeof(double);
    primitive_exponents = (double*) MemAllocator::allocate(block_size);
    contraction_coefficients = (double*) MemAllocator::allocate(block_size);
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
    ierr = xcint_set_functional(xcint_context, "b3lyp");
    ierr = xcint_set_basis(xcint_context,
                           XCINT_BASIS_SPHERICAL,
                           num_centers,
                           center_coordinates,
                           center_elements,
                           num_shells,
                           shell_centers,
                           shell_l_quantum_numbers,
                           shell_num_primitives,
                           primitive_exponents,
                           contraction_coefficients);

    int mat_dim = 19;
    double *dmat = NULL;
    double *xc_mat = NULL;
    block_size = mat_dim*mat_dim*sizeof(double);
    dmat = (double*) MemAllocator::allocate(block_size);
    std::fill(&dmat[0], &dmat[mat_dim*mat_dim], 0.0);
    xc_mat = (double*) MemAllocator::allocate(block_size);
    dmat[0] =   1.004012708363e+00;
    dmat[1] =   4.510701713634e-03;
    dmat[2] =  -1.144447782318e-02;
    dmat[3] =  -7.460192411188e-03;
    dmat[6] =   2.064755862872e-03;
    dmat[11] =  -9.392480166109e-05;
    dmat[13] =   1.626825285678e-04;
    dmat[14] =   3.997745549303e-03;
    dmat[15] =  -3.957074064980e-04;
    dmat[16] =  -4.898578168207e-04;
    dmat[19] =   4.510701713634e-03;
    dmat[20] =   9.194882886389e-01;
    dmat[21] =  -1.927744758881e-02;
    dmat[22] =   7.627553895467e-02;
    dmat[25] =   4.461232796069e-02;
    dmat[30] =  -1.407163850393e-03;
    dmat[32] =   2.437279283440e-03;
    dmat[33] =   1.516892605146e-01;
    dmat[34] =  -7.579356306024e-02;
    dmat[35] =   3.450787865374e-02;
    dmat[38] =  -1.144447782318e-02;
    dmat[39] =  -1.927744758881e-02;
    dmat[40] =   3.921259810966e-02;
    dmat[41] =   1.594594914555e-01;
    dmat[44] =  -2.180960509113e-02;
    dmat[49] =   1.737790248178e-03;
    dmat[51] =  -3.009941002742e-03;
    dmat[52] =  -1.480858033152e-01;
    dmat[53] =   6.494868921324e-02;
    dmat[54] =  -1.106166966204e-02;
    dmat[57] =  -7.460192411188e-03;
    dmat[58] =   7.627553895467e-02;
    dmat[59] =   1.594594914555e-01;
    dmat[60] =   6.762767650139e-01;
    dmat[63] =  -8.309659505673e-02;
    dmat[68] =   6.989016997025e-03;
    dmat[70] =  -1.210533253381e-02;
    dmat[71] =  -5.903237520676e-01;
    dmat[72] =   2.573872675260e-01;
    dmat[73] =  -4.018642354126e-02;
    dmat[80] =   9.296383420577e-01;
    dmat[83] =   2.687064170209e-02;
    dmat[85] =  -1.267372228454e-02;
    dmat[93] =   3.668493529941e-02;
    dmat[100] =   9.296383420577e-01;
    dmat[103] =   2.687064170209e-02;
    dmat[107] =  -1.267372228454e-02;
    dmat[113] =   3.668493529941e-02;
    dmat[114] =   2.064755862872e-03;
    dmat[115] =   4.461232796069e-02;
    dmat[116] =  -2.180960509113e-02;
    dmat[117] =  -8.309659505673e-02;
    dmat[120] =   1.341052289677e-02;
    dmat[125] =  -9.888841932932e-04;
    dmat[127] =   1.712797665585e-03;
    dmat[128] =   8.546811623303e-02;
    dmat[129] =  -3.783606332525e-02;
    dmat[130] =   7.250661042343e-03;
    dmat[137] =   2.687064170209e-02;
    dmat[140] =   7.766798687366e-04;
    dmat[142] =  -3.663263821346e-04;
    dmat[150] =   1.060356170457e-03;
    dmat[157] =   2.687064170209e-02;
    dmat[160] =   7.766798687367e-04;
    dmat[164] =  -3.663263821346e-04;
    dmat[170] =   1.060356170457e-03;
    dmat[175] =  -1.267372228454e-02;
    dmat[178] =  -3.663263821346e-04;
    dmat[180] =   1.727803483127e-04;
    dmat[188] =  -5.001242536768e-04;
    dmat[209] =  -9.392480166109e-05;
    dmat[210] =  -1.407163850393e-03;
    dmat[211] =   1.737790248178e-03;
    dmat[212] =   6.989016997025e-03;
    dmat[215] =  -9.888841932932e-04;
    dmat[220] =   7.751986552912e-05;
    dmat[222] =  -1.342683456924e-04;
    dmat[223] =  -6.626796710079e-03;
    dmat[224] =   2.912618333097e-03;
    dmat[225] =  -5.094002105433e-04;
    dmat[233] =  -1.267372228454e-02;
    dmat[236] =  -3.663263821346e-04;
    dmat[240] =   1.727803483127e-04;
    dmat[246] =  -5.001242536767e-04;
    dmat[247] =   1.626825285678e-04;
    dmat[248] =   2.437279283440e-03;
    dmat[249] =  -3.009941002742e-03;
    dmat[250] =  -1.210533253381e-02;
    dmat[253] =   1.712797665585e-03;
    dmat[258] =  -1.342683456924e-04;
    dmat[260] =   2.325595965874e-04;
    dmat[261] =   1.147794859329e-02;
    dmat[262] =  -5.044802935980e-03;
    dmat[263] =   8.823070460469e-04;
    dmat[266] =   3.997745549303e-03;
    dmat[267] =   1.516892605146e-01;
    dmat[268] =  -1.480858033152e-01;
    dmat[269] =  -5.903237520676e-01;
    dmat[272] =   8.546811623303e-02;
    dmat[277] =  -6.626796710079e-03;
    dmat[279] =   1.147794859329e-02;
    dmat[280] =   5.676124150509e-01;
    dmat[281] =  -2.498038873499e-01;
    dmat[282] =   4.443826475806e-02;
    dmat[285] =  -3.957074064980e-04;
    dmat[286] =  -7.579356306024e-02;
    dmat[287] =   6.494868921324e-02;
    dmat[288] =   2.573872675260e-01;
    dmat[291] =  -3.783606332525e-02;
    dmat[296] =   2.912618333097e-03;
    dmat[298] =  -5.044802935980e-03;
    dmat[299] =  -2.498038873499e-01;
    dmat[300] =   1.100325009820e-01;
    dmat[301] =  -1.979100163083e-02;
    dmat[304] =  -4.898578168207e-04;
    dmat[305] =   3.450787865374e-02;
    dmat[306] =  -1.106166966204e-02;
    dmat[307] =  -4.018642354126e-02;
    dmat[310] =   7.250661042343e-03;
    dmat[315] =  -5.094002105433e-04;
    dmat[317] =   8.823070460469e-04;
    dmat[318] =   4.443826475806e-02;
    dmat[319] =  -1.979100163083e-02;
    dmat[320] =   4.062611416622e-03;
    dmat[327] =   3.668493529941e-02;
    dmat[330] =   1.060356170457e-03;
    dmat[332] =  -5.001242536768e-04;
    dmat[340] =   1.447643042501e-03;
    dmat[347] =   3.668493529941e-02;
    dmat[350] =   1.060356170457e-03;
    dmat[354] =  -5.001242536767e-04;
    dmat[360] =   1.447643042501e-03;
    double xc_energy = 0.0;
    double num_electrons = 0.0;
    int dmat_to_pert[1]  = {0};
    int dmat_to_comp[1]  = {0};

    // generate grid
    numgrid_context_t *numgrid_context = numgrid_new();
    double radial_precision = 1.0e-12;
    int min_num_angular_points = 86;
    int max_num_angular_points = 302;
    int num_outer_centers = 0;
    double *outer_center_coordinates = NULL;
    int *outer_center_elements = NULL;
    ierr = numgrid_generate(numgrid_context,
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
    int num_points = numgrid_get_num_points(numgrid_context);
    double *grid_pw = (double*) numgrid_get_grid(numgrid_context);

    xcint_integrate(xcint_context,
                    XCINT_MODE_RKS,
                    num_points,
                    grid_pw,
                    0,
                    0,
                    0,
                    1,
                    dmat_to_pert,
                    dmat_to_comp,
                    dmat,
                    false,
                    xc_energy,
                    true,
                    xc_mat,
                    num_electrons);

    // free grid
    numgrid_free(numgrid_context);

    if (fabs(num_electrons - 9.999992182101e+00) > 1.0e-11) return_code++;
    double dot = 0.0;
    for (int i = 0; i < mat_dim*mat_dim; i++)
    {
        dot += xc_mat[i]*dmat[i];
    }
    if (fabs(dot - -5.620023209930e+00) > 1.0e-11) return_code++;
    MemAllocator::deallocate(dmat);
    MemAllocator::deallocate(xc_mat);
    MemAllocator::deallocate(center_coordinates);
    MemAllocator::deallocate(center_elements);
    MemAllocator::deallocate(shell_centers);
    MemAllocator::deallocate(shell_l_quantum_numbers);
    MemAllocator::deallocate(shell_num_primitives);
    MemAllocator::deallocate(primitive_exponents);
    MemAllocator::deallocate(contraction_coefficients);

    xcint_free(xcint_context);

    return return_code;
}
