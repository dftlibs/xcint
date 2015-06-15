#include <cmath>

#include "XCint.h"
#include "xcint.h"
#include "xcint_c_parameters.h"
#include "MemAllocator.h"
#include "numgrid.h"

int main(int argc, char** argv)
{
    int return_code = 0;

    xcint_context_t *xcint_context = xcint_new();

    int num_centers;
    int num_shells;
    double *center_xyz = NULL;
    int *center_element = NULL;
    int *shell_center = NULL;
    int *l_quantum_num = NULL;
    int *shell_num_primitives = NULL;
    double *primitive_exp = NULL;
    double *contraction_coef = NULL;
    size_t block_size;
    num_centers = 2;
    num_shells = 9;
    block_size = 3*num_centers*sizeof(double);
    center_xyz = (double*) MemAllocator::allocate(block_size);
    center_xyz[0] =   1.700000000000e+00;
    center_xyz[1] =   0.000000000000e+00;
    center_xyz[2] =   0.000000000000e+00;
    center_xyz[3] =   0.000000000000e+00;
    center_xyz[4] =   0.000000000000e+00;
    center_xyz[5] =   0.000000000000e+00;
    block_size = num_centers*sizeof(int);
    center_element = (int*) MemAllocator::allocate(block_size);
    center_element[0] = 9;
    center_element[1] = 1;
    block_size = num_shells*sizeof(int);
    shell_center = (int*) MemAllocator::allocate(block_size);
    l_quantum_num = (int*) MemAllocator::allocate(block_size);
    shell_num_primitives = (int*) MemAllocator::allocate(block_size);
    shell_center[0] = 1;
    l_quantum_num[0] = 0;
    shell_num_primitives[0] = 9;
    shell_center[1] = 1;
    l_quantum_num[1] = 0;
    shell_num_primitives[1] = 9;
    shell_center[2] = 1;
    l_quantum_num[2] = 0;
    shell_num_primitives[2] = 1;
    shell_center[3] = 1;
    l_quantum_num[3] = 1;
    shell_num_primitives[3] = 4;
    shell_center[4] = 1;
    l_quantum_num[4] = 1;
    shell_num_primitives[4] = 1;
    shell_center[5] = 1;
    l_quantum_num[5] = 2;
    shell_num_primitives[5] = 1;
    shell_center[6] = 2;
    l_quantum_num[6] = 0;
    shell_num_primitives[6] = 4;
    shell_center[7] = 2;
    l_quantum_num[7] = 0;
    shell_num_primitives[7] = 1;
    shell_center[8] = 2;
    l_quantum_num[8] = 1;
    shell_num_primitives[8] = 1;
    block_size = 31*sizeof(double);
    primitive_exp = (double*) MemAllocator::allocate(block_size);
    contraction_coef = (double*) MemAllocator::allocate(block_size);
    primitive_exp[0] =   1.471000000000e+04;
    contraction_coef[0] =   6.863650000000e-01;
    primitive_exp[1] =   2.207000000000e+03;
    contraction_coef[1] =   1.274350000000e+00;
    primitive_exp[2] =   5.028000000000e+02;
    contraction_coef[2] =   2.139130000000e+00;
    primitive_exp[3] =   1.426000000000e+02;
    contraction_coef[3] =   3.130550000000e+00;
    primitive_exp[4] =   4.647000000000e+01;
    contraction_coef[4] =   3.638230000000e+00;
    primitive_exp[5] =   1.670000000000e+01;
    contraction_coef[5] =   2.641480000000e+00;
    primitive_exp[6] =   6.356000000000e+00;
    contraction_coef[6] =   7.553570000000e-01;
    primitive_exp[7] =   1.316000000000e+00;
    contraction_coef[7] =   1.342700000000e-02;
    primitive_exp[8] =   3.897000000000e-01;
    contraction_coef[8] =  -8.197600000000e-04;
    primitive_exp[9] =   1.471000000000e+04;
    contraction_coef[9] =  -1.570740000000e-01;
    primitive_exp[10] =   2.207000000000e+03;
    contraction_coef[10] =  -3.001720000000e-01;
    primitive_exp[11] =   5.028000000000e+02;
    contraction_coef[11] =  -4.915140000000e-01;
    primitive_exp[12] =   1.426000000000e+02;
    contraction_coef[12] =  -7.849910000000e-01;
    primitive_exp[13] =   4.647000000000e+01;
    contraction_coef[13] =  -9.347560000000e-01;
    primitive_exp[14] =   1.670000000000e+01;
    contraction_coef[14] =  -1.005480000000e+00;
    primitive_exp[15] =   6.356000000000e+00;
    contraction_coef[15] =  -3.204660000000e-01;
    primitive_exp[16] =   1.316000000000e+00;
    contraction_coef[16] =   4.928530000000e-01;
    primitive_exp[17] =   3.897000000000e-01;
    contraction_coef[17] =   1.999410000000e-01;
    primitive_exp[18] =   3.897000000000e-01;
    contraction_coef[18] =   3.515260000000e-01;
    primitive_exp[19] =   2.267000000000e+01;
    contraction_coef[19] =   3.164380000000e+00;
    primitive_exp[20] =   4.977000000000e+00;
    contraction_coef[20] =   2.497710000000e+00;
    primitive_exp[21] =   1.347000000000e+00;
    contraction_coef[21] =   1.051860000000e+00;
    primitive_exp[22] =   3.471000000000e-01;
    contraction_coef[22] =   1.739750000000e-01;
    primitive_exp[23] =   3.471000000000e-01;
    contraction_coef[23] =   3.797590000000e-01;
    primitive_exp[24] =   1.640000000000e+00;
    contraction_coef[24] =   6.775590000000e+00;
    primitive_exp[25] =   1.301000000000e+01;
    contraction_coef[25] =   9.610660000000e-02;
    primitive_exp[26] =   1.962000000000e+00;
    contraction_coef[26] =   1.630200000000e-01;
    primitive_exp[27] =   4.446000000000e-01;
    contraction_coef[27] =   1.855450000000e-01;
    primitive_exp[28] =   1.220000000000e-01;
    contraction_coef[28] =   7.374380000000e-02;
    primitive_exp[29] =   1.220000000000e-01;
    contraction_coef[29] =   1.471230000000e-01;
    primitive_exp[30] =   7.270000000000e-01;
    contraction_coef[30] =   9.568810000000e-01;
    double hfx, mu, beta; // we don't care about it here
    xcint_set_functional(xcint_context, "lda", hfx, mu, beta);
    xcint_set_basis(xcint_context,
                    XCINT_BASIS_SPHERICAL,
                    num_centers,
                    center_xyz,
                    center_element,
                    num_shells,
                    shell_center,
                    l_quantum_num,
                    shell_num_primitives,
                    primitive_exp,
                    contraction_coef);

    int mat_dim = 19;
    double *dmat = NULL;
    double *xc_mat = NULL;
    block_size = mat_dim*mat_dim*sizeof(double);
    dmat = (double*) MemAllocator::allocate(block_size);
    std::fill(&dmat[0], &dmat[mat_dim*mat_dim], 0.0);
    xc_mat = (double*) MemAllocator::allocate(block_size);
    dmat[0] =   1.007006744703e+00;
    dmat[1] =   8.524020246319e-03;
    dmat[2] =  -1.913517463803e-02;
    dmat[3] =  -8.640040075696e-03;
    dmat[6] =   2.016135048780e-03;
    dmat[11] =  -1.373086358909e-04;
    dmat[13] =   2.378255336812e-04;
    dmat[14] =   3.313510584148e-03;
    dmat[15] =   3.582318408589e-04;
    dmat[16] =  -9.808679472161e-04;
    dmat[19] =   8.524020246319e-03;
    dmat[20] =   8.864234480930e-01;
    dmat[21] =  -4.073012026704e-03;
    dmat[22] =   7.024766747215e-02;
    dmat[25] =   4.742225668967e-02;
    dmat[30] =  -1.820026289381e-03;
    dmat[32] =   3.152378004318e-03;
    dmat[33] =   1.513601817934e-01;
    dmat[34] =  -7.276087881061e-02;
    dmat[35] =   3.476210852707e-02;
    dmat[38] =  -1.913517463803e-02;
    dmat[39] =  -4.073012026704e-03;
    dmat[40] =   3.808645776634e-02;
    dmat[41] =   1.599421660580e-01;
    dmat[44] =  -2.061395923696e-02;
    dmat[49] =   1.679385562496e-03;
    dmat[51] =  -2.908781119740e-03;
    dmat[52] =  -1.408831586336e-01;
    dmat[53] =   6.131376191995e-02;
    dmat[54] =  -9.198561533752e-03;
    dmat[57] =  -8.640040075696e-03;
    dmat[58] =   7.024766747215e-02;
    dmat[59] =   1.599421660580e-01;
    dmat[60] =   6.853469832411e-01;
    dmat[63] =  -8.272711349281e-02;
    dmat[68] =   6.941961531658e-03;
    dmat[70] =  -1.202383007702e-02;
    dmat[71] =  -5.830734532642e-01;
    dmat[72] =   2.532122504427e-01;
    dmat[73] =  -3.571545507150e-02;
    dmat[80] =   9.225234222380e-01;
    dmat[83] =   3.022090769501e-02;
    dmat[85] =  -1.311210715170e-02;
    dmat[93] =   3.850540048750e-02;
    dmat[100] =   9.225234222380e-01;
    dmat[103] =   3.022090769502e-02;
    dmat[107] =  -1.311210715170e-02;
    dmat[113] =   3.850540048750e-02;
    dmat[114] =   2.016135048780e-03;
    dmat[115] =   4.742225668967e-02;
    dmat[116] =  -2.061395923696e-02;
    dmat[117] =  -8.272711349281e-02;
    dmat[120] =   1.354027328425e-02;
    dmat[125] =  -9.989208178770e-04;
    dmat[127] =   1.730181609301e-03;
    dmat[128] =   8.380299335376e-02;
    dmat[129] =  -3.683950216674e-02;
    dmat[130] =   6.753429473392e-03;
    dmat[137] =   3.022090769501e-02;
    dmat[140] =   9.900055000170e-04;
    dmat[142] =  -4.295389909529e-04;
    dmat[150] =   1.261396866292e-03;
    dmat[157] =   3.022090769502e-02;
    dmat[160] =   9.900055000176e-04;
    dmat[164] =  -4.295389909531e-04;
    dmat[170] =   1.261396866293e-03;
    dmat[175] =  -1.311210715170e-02;
    dmat[178] =  -4.295389909529e-04;
    dmat[180] =   1.863663835663e-04;
    dmat[188] =  -5.472890172115e-04;
    dmat[209] =  -1.373086358909e-04;
    dmat[210] =  -1.820026289381e-03;
    dmat[211] =   1.679385562496e-03;
    dmat[212] =   6.941961531658e-03;
    dmat[215] =  -9.989208178770e-04;
    dmat[220] =   7.760570466053e-05;
    dmat[222] =  -1.344170234292e-04;
    dmat[223] =  -6.513776896370e-03;
    dmat[224] =   2.848946299205e-03;
    dmat[225] =  -4.723568425938e-04;
    dmat[233] =  -1.311210715170e-02;
    dmat[236] =  -4.295389909531e-04;
    dmat[240] =   1.863663835662e-04;
    dmat[246] =  -5.472890172115e-04;
    dmat[247] =   2.378255336812e-04;
    dmat[248] =   3.152378004318e-03;
    dmat[249] =  -2.908781119740e-03;
    dmat[250] =  -1.202383007702e-02;
    dmat[253] =   1.730181609301e-03;
    dmat[258] =  -1.344170234292e-04;
    dmat[260] =   2.328171139816e-04;
    dmat[261] =   1.128219253368e-02;
    dmat[262] =  -4.934519738259e-03;
    dmat[263] =   8.181460506753e-04;
    dmat[266] =   3.313510584148e-03;
    dmat[267] =   1.513601817934e-01;
    dmat[268] =  -1.408831586336e-01;
    dmat[269] =  -5.830734532642e-01;
    dmat[272] =   8.380299335376e-02;
    dmat[277] =  -6.513776896370e-03;
    dmat[279] =   1.128219253368e-02;
    dmat[280] =   5.467981651861e-01;
    dmat[281] =  -2.391585905977e-01;
    dmat[282] =   3.962356807742e-02;
    dmat[285] =   3.582318408589e-04;
    dmat[286] =  -7.276087881061e-02;
    dmat[287] =   6.131376191995e-02;
    dmat[288] =   2.532122504427e-01;
    dmat[291] =  -3.683950216674e-02;
    dmat[296] =   2.848946299205e-03;
    dmat[298] =  -4.934519738259e-03;
    dmat[299] =  -2.391585905977e-01;
    dmat[300] =   1.046575883118e-01;
    dmat[301] =  -1.751801112317e-02;
    dmat[304] =  -9.808679472161e-04;
    dmat[305] =   3.476210852707e-02;
    dmat[306] =  -9.198561533752e-03;
    dmat[307] =  -3.571545507150e-02;
    dmat[310] =   6.753429473392e-03;
    dmat[315] =  -4.723568425938e-04;
    dmat[317] =   8.181460506753e-04;
    dmat[318] =   3.962356807742e-02;
    dmat[319] =  -1.751801112317e-02;
    dmat[320] =   3.543728922920e-03;
    dmat[327] =   3.850540048750e-02;
    dmat[330] =   1.261396866292e-03;
    dmat[332] =  -5.472890172115e-04;
    dmat[340] =   1.607185065400e-03;
    dmat[347] =   3.850540048750e-02;
    dmat[350] =   1.261396866293e-03;
    dmat[354] =  -5.472890172115e-04;
    dmat[360] =   1.607185065400e-03;
    double xc_energy = 0.0;
    double num_electrons = 0.0;
    int dmat_to_pert[1]  = {0};
    int dmat_to_comp[1]  = {0};

    // generate grid
    numgrid_context_t *numgrid_context = numgrid_new();
    double radial_precision = 1.0e-12;
    int angular_min = 86;
    int angular_max = 302;
    int num_outer_centers = 0;
    double *outer_center_xyz = NULL;
    int *outer_center_element = NULL;
    numgrid_generate(numgrid_context,
                     radial_precision,
                     angular_min,
                     angular_max,
                     num_centers,
                     center_xyz,
                     center_element,
                     num_outer_centers,
                     outer_center_xyz,
                     outer_center_element,
                     num_shells,
                     shell_center,
                     l_quantum_num,
                     shell_num_primitives,
                     primitive_exp);
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

    if (fabs(num_electrons - 9.999992074832e+00) > 1.0e-11) return_code++;
    double dot = 0.0;
    for (int i = 0; i < mat_dim*mat_dim; i++)
    {
        dot += xc_mat[i]*dmat[i];
    }
    if (fabs(dot - -6.729996811122e+00) > 1.0e-11) return_code++;
    MemAllocator::deallocate(dmat);
    MemAllocator::deallocate(xc_mat);
    MemAllocator::deallocate(center_xyz);
    MemAllocator::deallocate(center_element);
    MemAllocator::deallocate(shell_center);
    MemAllocator::deallocate(l_quantum_num);
    MemAllocator::deallocate(shell_num_primitives);
    MemAllocator::deallocate(primitive_exp);
    MemAllocator::deallocate(contraction_coef);

    xcint_free(xcint_context);

    return return_code;
}
