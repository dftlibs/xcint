#include "integrator.h"

#include "compress.h"
#include "density.h"
#include "generated_parameters.h"
#include "xcint_parameters.h"

#include "xcfun.h"
#include "xcint.h"

#include <cassert>
#include <cmath>
#include <ctime>

#include <algorithm>
#include <cstdlib>
#include <fstream>

#define AS_TYPE(Type, Obj) reinterpret_cast<Type *>(Obj)
#define AS_CTYPE(Type, Obj) reinterpret_cast<const Type *>(Obj)

// FIXME ugly global
xc_functional xcfun;
int dens_offset;

XCINT_API
xcint_context_t *xcint_new_context()
{
    return AS_TYPE(xcint_context_t, new XCint());
}
XCint::XCint() { balboa_context = balboa_new_context(); }

XCINT_API
void xcint_free_context(xcint_context_t *xcint_context)
{
    if (!xcint_context)
        return;
    delete AS_TYPE(XCint, xcint_context);
}
XCint::~XCint() { balboa_free_context(balboa_context); }

XCINT_API
int xcint_set_functional(xcint_context_t *context, const char *line)
{
    return AS_TYPE(XCint, context)->set_functional(line);
}
int XCint::set_functional(const char *line)
{
    fun.set_functional(line);
    return 0;
}

XCINT_API
int xcint_set_basis(xcint_context_t *context,
                    const xcint_basis_t basis_type,
                    const int num_centers,
                    const double center_coordinates[],
                    const int num_shells,
                    const int shell_centers[],
                    const int shell_l_quantum_numbers[],
                    const int shell_num_primitives[],
                    const double primitive_exponents[],
                    const double contraction_coefficients[])
{
    return AS_TYPE(XCint, context)
        ->set_basis(basis_type,
                    num_centers,
                    center_coordinates,
                    num_shells,
                    shell_centers,
                    shell_l_quantum_numbers,
                    shell_num_primitives,
                    primitive_exponents,
                    contraction_coefficients);
}
int XCint::set_basis(const int basis_type,
                     const int num_centers,
                     const double center_coordinates_bohr[],
                     const int num_shells,
                     const int shell_centers[],
                     const int shell_l_quantum_numbers[],
                     const int shell_num_primitives[],
                     const double primitive_exponents[],
                     const double contraction_coefficients[])
{
    int ierr = balboa_set_basis(balboa_context,
                                basis_type,
                                num_centers,
                                center_coordinates_bohr,
                                num_shells,
                                shell_centers,
                                shell_l_quantum_numbers,
                                shell_num_primitives,
                                primitive_exponents,
                                contraction_coefficients);

    return ierr;
}

void XCint::integrate_batch(const double dmat[],
                            const bool get_exc,
                            double &exc,
                            const bool get_vxc,
                            double vxc[],
                            double &num_electrons,
                            const int geo_coor[],
                            const bool use_dmat[],
                            const int ipoint,
                            const int geo_derv_order,
                            const int max_ao_order_g,
                            const int block_length,
                            const int num_variables,
                            const int num_perturbations,
                            const int num_fields,
                            const int mat_dim,
                            const bool get_gradient,
                            const bool get_tau,
                            const int dmat_index[],
                            const double grid_x_bohr[],
                            const double grid_y_bohr[],
                            const double grid_z_bohr[],
                            const double grid_w[])
//  const double grid_w[]) const
{
    double *n = new double[AO_BLOCK_LENGTH * num_variables * MAX_NUM_DENSITIES];
    double *u = new double[AO_BLOCK_LENGTH * num_variables];

    std::vector<int> coor;
    double prefactors[5] = {1.0, 2.0, 2.0, 2.0, 0.5};

    bool n_is_used[MAX_NUM_DENSITIES];
    std::fill(&n_is_used[0], &n_is_used[MAX_NUM_DENSITIES], false);

    int max_ao_geo_order = max_ao_order_g; // FIXME

    int buffer_len = balboa_get_buffer_len(
        balboa_context, max_ao_geo_order, AO_BLOCK_LENGTH);
    //  FIXME should be:
    //  int buffer_len = balboa_get_buffer_len(balboa_context, max_ao_geo_order,
    //  block_length);

    ao_length = buffer_len;
    ao = new double[buffer_len];

    std::fill(&ao[0], &ao[buffer_len], 0.0);

    int ierr = balboa_get_ao(balboa_context,
                             max_ao_geo_order,
                             block_length,
                             &grid_x_bohr[ipoint],
                             &grid_y_bohr[ipoint],
                             &grid_z_bohr[ipoint],
                             ao);

    if (!n_is_used[0])
    {
        std::fill(&n[0], &n[block_length * num_variables], 0.0);
        n_is_used[0] = true;
    }

    double *ao_compressed = new double[buffer_len];
    int *ao_compressed_index = new int[buffer_len];
    int ao_compressed_num;
    int num_aos = balboa_get_num_aos(balboa_context);
    int *ao_centers = new int[num_aos];
    for (int i = 0; i < num_aos; i++)
    {
        ao_centers[i] = balboa_get_ao_center(balboa_context, i);
    }
    int slice_offsets[4];
    auto get_geo_offset = [&](int i, int j, int k) {
        return balboa_get_geo_offset(balboa_context, i, j, k);
    };

    compute_slice_offsets(std::vector<int>(), slice_offsets);
    compress(get_gradient,
             block_length,
             ao_compressed_num,
             ao_compressed_index,
             ao_compressed,
             num_aos,
             ao,
             ao_centers,
             std::vector<int>(),
             slice_offsets);
    get_density(mat_dim,
                block_length,
                get_gradient,
                get_tau,
                prefactors,
                n,
                dmat,
                true,
                true,
                ao_compressed_num,
                ao_compressed_index,
                ao_compressed,
                ao_compressed_num,
                ao_compressed_index,
                ao_compressed);

    for (int ib = 0; ib < block_length; ib++)
    {
        num_electrons += grid_w[ipoint + ib] * n[ib];
    }

    // expectation value contribution
    if (get_exc)
    {

        double *xcin = NULL;
        double *xcout = NULL;

        assert(geo_derv_order < 5);

        //                       *
        //                       *  *
        //                       *  *  *  *
        //                       *  *  *  *  *  *  *  *
        //                       *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
        //           i j k l     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
        //   ------------------------------------------------------------------
        //    0 0                0
        //    1 i    *           1  0
        //    2 j      *         2     0
        //    3 ij   * *         3  2  1  0
        //    4 k        *       4           0
        //    5 ik   *   *       5  4        1  0
        //    6 jk     * *       6     4     2     0
        //    7 ijk  * * *       7  6  5  4  3  2  1  0
        //    8 l          *     8                       0
        //    9 il   *     *     9  8                    1  0
        //   10 jl     *   *    10     8                 2     0
        //   11 ijl  * *   *    11 10  9  8              3  2  1  0
        //   12 kl       * *    12           8           4           0
        //   13 ikl  *   * *    13 12        9  8        5  4        1  0
        //   14 jkl    * * *    14    12    10     8     6     4     2     0
        //   15 ijkl * * * *    15 14 13 12 11 10  9  8  7  6  5  4  3  2  1  0

        // this takes care of the first column
        for (int k = 1; k < MAX_NUM_DENSITIES; k++)
        {
            if (use_dmat[k])
            {
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                get_density(
                    mat_dim,
                    block_length,
                    get_gradient,
                    get_tau,
                    prefactors,
                    &n[k * block_length * num_variables],
                    &dmat[dmat_index[k]],
                    false,
                    false, // FIXME can be true based on dmat, saving possible
                    ao_compressed_num,
                    ao_compressed_index,
                    ao_compressed,
                    ao_compressed_num,
                    ao_compressed_index,
                    ao_compressed);
            }
        }

#include "ave_contributions.h"

        dens_offset = fun.set_order(num_perturbations, xcfun);

        xcin = new double[num_variables * dens_offset * block_length];
        std::fill(
            &xcin[0], &xcin[num_variables * dens_offset * block_length], 0.0);

        xcout = new double[dens_offset * block_length];

        for (int k = 0; k < MAX_NUM_DENSITIES; k++)
        {
            if (n_is_used[k])
            {
                for (int ivar = 0; ivar < num_variables; ivar++)
                {
                    for (int ib = 0; ib < block_length; ib++)
                    {
                        xcin[ib * num_variables * dens_offset +
                             ivar * dens_offset + k] =
                            n[k * block_length * num_variables +
                              ivar * block_length + ib];
                    }
                }
            }
        }

        double sum = 0.0;
        for (int ib = 0; ib < block_length; ib++)
        {
            if (n[ib] > 1.0e-14 and std::abs(grid_w[ipoint + ib]) > 1.0e-30)
            {
                xc_eval(xcfun,
                        &xcin[ib * num_variables * dens_offset],
                        &xcout[ib * dens_offset]);
                sum += xcout[ib * dens_offset + dens_offset - 1] *
                       grid_w[ipoint + ib];
            }
        }
        exc += sum;

        delete[] xcin;
        delete[] xcout;
    }

    // matrix contribution
    if (get_vxc)
    {
        int k;

        bool contribution_is_implemented = false;

        if (geo_derv_order == 0) // no geo dervs
        {
            contribution_is_implemented = true;

            for (int ifield = 0; ifield < num_fields; ifield++)
            {
                k = (int)pow(2, ifield); // FIXME double check this
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                get_density(
                    mat_dim,
                    block_length,
                    get_gradient,
                    get_tau,
                    prefactors,
                    &n[k * block_length * num_variables],
                    &dmat[(ifield + 1) * mat_dim * mat_dim],
                    false,
                    false, // FIXME can be true depending on perturbation
                           // (savings possible)
                    ao_compressed_num,
                    ao_compressed_index,
                    ao_compressed,
                    ao_compressed_num,
                    ao_compressed_index,
                    ao_compressed);
            }

            distribute_matrix2(block_length,
                               num_variables,
                               num_perturbations,
                               mat_dim,
                               prefactors,
                               ipoint,
                               n_is_used,
                               n,
                               u,
                               vxc,
                               exc,
                               coor,
                               grid_w);
        }

        if (geo_derv_order > 0) // we have geo dervs
        {
            if (geo_derv_order == 1 && num_fields == 0)
            {
                contribution_is_implemented = true;

                // M    d_nn  n_i
                k = 1;
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                get_dens_geo_derv(mat_dim,
                                  num_aos,
                                  block_length,
                                  buffer_len,
                                  ao,
                                  ao_centers,
                                  get_gradient,
                                  get_tau,
                                  coor,
                                  get_geo_offset,
                                  &n[k * block_length * num_variables],
                                  &dmat[0]);
                coor.clear();
                distribute_matrix2(block_length,
                                   num_variables,
                                   1,
                                   mat_dim,
                                   prefactors,
                                   ipoint,
                                   n_is_used,
                                   n,
                                   u,
                                   vxc,
                                   exc,
                                   coor,
                                   grid_w);
                n_is_used[1] = false;

                // M_i  d_n
                coor.push_back(geo_coor[0]);
                distribute_matrix2(block_length,
                                   num_variables,
                                   0,
                                   mat_dim,
                                   prefactors,
                                   ipoint,
                                   n_is_used,
                                   n,
                                   u,
                                   vxc,
                                   exc,
                                   coor,
                                   grid_w);
                coor.clear();
            }

            if (geo_derv_order == 2 && num_fields == 0)
            {
                contribution_is_implemented = true;

                // M_ij d_n
                coor.push_back(geo_coor[0]);
                coor.push_back(geo_coor[1]);
                distribute_matrix2(block_length,
                                   num_variables,
                                   0,
                                   mat_dim,
                                   prefactors,
                                   ipoint,
                                   n_is_used,
                                   n,
                                   u,
                                   vxc,
                                   exc,
                                   coor,
                                   grid_w);
                coor.clear();

                // FIXME add shortcut if i == j
                // M_i  d_nn  n_j
                k = 1;
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[1]);
                get_dens_geo_derv(mat_dim,
                                  num_aos,
                                  block_length,
                                  buffer_len,
                                  ao,
                                  ao_centers,
                                  get_gradient,
                                  get_tau,
                                  coor,
                                  get_geo_offset,
                                  &n[k * block_length * num_variables],
                                  &dmat[0]);
                coor.clear();
                coor.push_back(geo_coor[0]);
                distribute_matrix2(block_length,
                                   num_variables,
                                   1,
                                   mat_dim,
                                   prefactors,
                                   ipoint,
                                   n_is_used,
                                   n,
                                   u,
                                   vxc,
                                   exc,
                                   coor,
                                   grid_w);
                coor.clear();
                n_is_used[1] = false;

                // M_j  d_nn  n_i
                k = 1;
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                get_dens_geo_derv(mat_dim,
                                  num_aos,
                                  block_length,
                                  buffer_len,
                                  ao,
                                  ao_centers,
                                  get_gradient,
                                  get_tau,
                                  coor,
                                  get_geo_offset,
                                  &n[k * block_length * num_variables],
                                  &dmat[0]);
                coor.clear();
                coor.push_back(geo_coor[1]);
                distribute_matrix2(block_length,
                                   num_variables,
                                   1,
                                   mat_dim,
                                   prefactors,
                                   ipoint,
                                   n_is_used,
                                   n,
                                   u,
                                   vxc,
                                   exc,
                                   coor,
                                   grid_w);
                coor.clear();
                n_is_used[1] = false;

                // M    d_nnn n_i n_j
                // M    d_nn  n_ij
                // FIXME stupidly recalculating n_i and n_j
                k = 1;
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                get_dens_geo_derv(mat_dim,
                                  num_aos,
                                  block_length,
                                  buffer_len,
                                  ao,
                                  ao_centers,
                                  get_gradient,
                                  get_tau,
                                  coor,
                                  get_geo_offset,
                                  &n[k * block_length * num_variables],
                                  &dmat[0]);
                coor.clear();
                k = 2;
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[1]);
                get_dens_geo_derv(mat_dim,
                                  num_aos,
                                  block_length,
                                  buffer_len,
                                  ao,
                                  ao_centers,
                                  get_gradient,
                                  get_tau,
                                  coor,
                                  get_geo_offset,
                                  &n[k * block_length * num_variables],
                                  &dmat[0]);
                coor.clear();
                k = 3;
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                coor.push_back(geo_coor[1]);
                get_dens_geo_derv(mat_dim,
                                  num_aos,
                                  block_length,
                                  buffer_len,
                                  ao,
                                  ao_centers,
                                  get_gradient,
                                  get_tau,
                                  coor,
                                  get_geo_offset,
                                  &n[k * block_length * num_variables],
                                  &dmat[0]);
                coor.clear();
                distribute_matrix2(block_length,
                                   num_variables,
                                   2,
                                   mat_dim,
                                   prefactors,
                                   ipoint,
                                   n_is_used,
                                   n,
                                   u,
                                   vxc,
                                   exc,
                                   coor,
                                   grid_w);
                n_is_used[1] = false;
                n_is_used[2] = false;
                n_is_used[3] = false;
            }

            if (geo_derv_order == 1 && num_fields == 1)
            {
                contribution_is_implemented = true;

                // M_i  d_nn  n_a
                k = 1;
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                get_density(mat_dim,
                            block_length,
                            get_gradient,
                            get_tau,
                            prefactors,
                            &n[k * block_length * num_variables],
                            &dmat[1 * mat_dim * mat_dim],
                            false,
                            false, // FIXME can be true depending
                                   // on perturbation (savings
                                   // possible)
                            ao_compressed_num,
                            ao_compressed_index,
                            ao_compressed,
                            ao_compressed_num,
                            ao_compressed_index,
                            ao_compressed);
                coor.push_back(geo_coor[0]);
                distribute_matrix2(block_length,
                                   num_variables,
                                   1,
                                   mat_dim,
                                   prefactors,
                                   ipoint,
                                   n_is_used,
                                   n,
                                   u,
                                   vxc,
                                   exc,
                                   coor,
                                   grid_w);
                coor.clear();

                // M    d_nnn n_i n_a
                // M    d_nn  n_ia
                k = 2;
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                get_dens_geo_derv(mat_dim,
                                  num_aos,
                                  block_length,
                                  buffer_len,
                                  ao,
                                  ao_centers,
                                  get_gradient,
                                  get_tau,
                                  coor,
                                  get_geo_offset,
                                  &n[k * block_length * num_variables],
                                  &dmat[0]);
                coor.clear();
                k = 3;
                if (!n_is_used[k])
                {
                    std::fill(&n[k * block_length * num_variables],
                              &n[(k + 1) * block_length * num_variables],
                              0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                get_dens_geo_derv(mat_dim,
                                  num_aos,
                                  block_length,
                                  buffer_len,
                                  ao,
                                  ao_centers,
                                  get_gradient,
                                  get_tau,
                                  coor,
                                  get_geo_offset,
                                  &n[k * block_length * num_variables],
                                  &dmat[1 * mat_dim * mat_dim]);
                coor.clear();
                distribute_matrix2(block_length,
                                   num_variables,
                                   2,
                                   mat_dim,
                                   prefactors,
                                   ipoint,
                                   n_is_used,
                                   n,
                                   u,
                                   vxc,
                                   exc,
                                   coor,
                                   grid_w);
                n_is_used[1] = false;
                n_is_used[2] = false;
                n_is_used[3] = false;
            }
        }

        if (!contribution_is_implemented)
        {
            fprintf(stderr,
                    "ERROR: XCint matrix contribution for "
                    "geo_derv_order=%i and num_fields=%i not "
                    "implemented\n",
                    geo_derv_order,
                    num_fields);
            exit(-1);
        }
    }

    delete[] n;
    delete[] u;
    delete[] ao;
    delete[] ao_compressed;
    delete[] ao_compressed_index;
    delete[] ao_centers;
}

XCINT_API
// int xcint_integrate_scf(const xcint_context_t *context,
int xcint_integrate_scf(xcint_context_t *context,
                        const xcint_mode_t mode,
                        const int num_points,
                        const double grid_x_bohr[],
                        const double grid_y_bohr[],
                        const double grid_z_bohr[],
                        const double grid_w[],
                        const double dmat[],
                        double *exc,
                        double vxc[],
                        double *num_electrons)
{
    int num_perturbations = 0;
    xcint_perturbation_t *perturbations = NULL;
    int *components = NULL;
    int num_dmat = 1;
    int *perturbation_indices = NULL;
    bool get_exc = true;
    bool get_vxc = true;
    return AS_TYPE(XCint, context)
        ->integrate(mode,
                    num_points,
                    grid_x_bohr,
                    grid_y_bohr,
                    grid_z_bohr,
                    grid_w,
                    num_perturbations,
                    perturbations,
                    components,
                    num_dmat,
                    perturbation_indices,
                    dmat,
                    get_exc,
                    exc,
                    get_vxc,
                    vxc,
                    num_electrons);
}

XCINT_API
// int xcint_integrate(const xcint_context_t *context,
int xcint_integrate(xcint_context_t *context,
                    const xcint_mode_t mode,
                    const int num_points,
                    const double grid_x_bohr[],
                    const double grid_y_bohr[],
                    const double grid_z_bohr[],
                    const double grid_w[],
                    const int num_perturbations,
                    const xcint_perturbation_t perturbations[],
                    const int components[],
                    const int num_dmat,
                    const int perturbation_indices[],
                    const double dmat[],
                    const bool get_exc,
                    double *exc,
                    const bool get_vxc,
                    double vxc[],
                    double *num_electrons)
{
    return AS_TYPE(XCint, context)
        ->integrate(mode,
                    num_points,
                    grid_x_bohr,
                    grid_y_bohr,
                    grid_z_bohr,
                    grid_w,
                    num_perturbations,
                    perturbations,
                    components,
                    num_dmat,
                    perturbation_indices,
                    dmat,
                    get_exc,
                    exc,
                    get_vxc,
                    vxc,
                    num_electrons);
}
int XCint::integrate(const xcint_mode_t mode,
                     const int num_points,
                     const double grid_x_bohr[],
                     const double grid_y_bohr[],
                     const double grid_z_bohr[],
                     const double grid_w[],
                     const int num_perturbations,
                     const xcint_perturbation_t perturbations[],
                     const int components[],
                     const int num_dmat,
                     const int perturbation_indices[],
                     const double dmat[],
                     const bool get_exc,
                     double *exc,
                     const bool get_vxc,
                     double vxc[],
                     //  double *num_electrons) const
                     double *num_electrons)
{
    xcfun = xc_new_functional();
    for (size_t i = 0; i < fun.keys.size(); i++)
    {
        int ierr = xc_set(xcfun, fun.keys[i].c_str(), fun.weights[i]);
        if (ierr != 0)
        {
            fprintf(stderr,
                    "ERROR in fun init: \"%s\" not recognized, quitting.\n",
                    fun.keys[i].c_str());
            exit(-1);
        }
    }

    int num_variables;
    int max_ao_order_g;

    double a;

    std::vector<int> coor;

    int mat_dim = balboa_get_num_aos(balboa_context);

    int geo_derv_order = 0;
    int num_fields = 0;
    for (int i = 0; i < num_perturbations; i++)
    {
        if (perturbations[i] == XCINT_PERT_GEO)
            geo_derv_order++;
        if (perturbations[i] == XCINT_PERT_EL)
            num_fields++;
    }

    int *geo_coor = NULL;
    if (geo_derv_order > 0)
    {
        geo_coor = new int[geo_derv_order];
        for (int i = 0; i < num_perturbations; i++)
        {
            if (perturbations[i] == XCINT_PERT_GEO)
                geo_coor[i] = components[2 * i]; // FIXME
        }
    }

    assert(mode == XCINT_MODE_RKS);

    *exc = 0.0;
    *num_electrons = 0.0;

    if (get_vxc)
        std::fill(&vxc[0], &vxc[mat_dim * mat_dim], 0.0);

    bool get_gradient;
    bool get_tau;

    max_ao_order_g = 0;

    if (fun.is_tau_mgga)
    {
        num_variables = 5;
        get_gradient = true;
        get_tau = true;
        max_ao_order_g++;
    }
    else if (fun.is_gga)
    {
        num_variables = 4;
        get_gradient = true;
        get_tau = false;
        max_ao_order_g++;
    }
    else
    {
        num_variables = 1;
        get_gradient = false;
        get_tau = false;
    }

    for (int i = 0; i < geo_derv_order; i++)
    {
        max_ao_order_g++;
    }

    int *dmat_index = new int[MAX_NUM_DENSITIES];
    std::fill(&dmat_index[0], &dmat_index[MAX_NUM_DENSITIES], 0);

    bool *use_dmat = new bool[MAX_NUM_DENSITIES];
    std::fill(&use_dmat[0], &use_dmat[MAX_NUM_DENSITIES], false);

    assert(num_dmat <= MAX_NUM_DENSITIES);

    use_dmat[0] = true;
    dmat_index[0] = 0;
    for (int k = 1; k < num_dmat; k++)
    {
        use_dmat[perturbation_indices[k]] = true;
        dmat_index[perturbation_indices[k]] = k * mat_dim * mat_dim;
    }

    assert(num_perturbations < 7);

    int block_length;
    int num_points_left = num_points;
    int num_batches = num_points / AO_BLOCK_LENGTH;
    if (num_points % AO_BLOCK_LENGTH != 0)
        num_batches++;

    double exc_local = *exc;
    double num_electrons_local = *num_electrons;
    double *vxc_local = NULL;

    if (get_vxc)
        vxc_local = &vxc[0];

    for (int ibatch = 0; ibatch < num_batches; ibatch++)
    {
        int ipoint = ibatch * AO_BLOCK_LENGTH;

        if (num_points_left < AO_BLOCK_LENGTH)
        {
            block_length = num_points_left;
        }
        else
        {
            block_length = AO_BLOCK_LENGTH;
        }

        integrate_batch(dmat,
                        get_exc,
                        exc_local,
                        get_vxc,
                        vxc_local,
                        num_electrons_local,
                        geo_coor,
                        use_dmat,
                        ipoint,
                        geo_derv_order,
                        max_ao_order_g,
                        block_length,
                        num_variables,
                        num_perturbations,
                        num_fields,
                        mat_dim,
                        get_gradient,
                        get_tau,
                        dmat_index,
                        grid_x_bohr,
                        grid_y_bohr,
                        grid_z_bohr,
                        grid_w);

        num_points_left -= block_length;
    }

    *exc = exc_local;
    *num_electrons = num_electrons_local;

    delete[] use_dmat;
    delete[] dmat_index;
    delete[] geo_coor;

    if (get_vxc)
    {
        // symmetrize result matrix
        for (int k = 0; k < mat_dim; k++)
        {
            for (int l = 0; l < k; l++)
            {
                a = vxc[k * mat_dim + l] + vxc[l * mat_dim + k];
                vxc[k * mat_dim + l] = 0.5 * a;
                vxc[l * mat_dim + k] = 0.5 * a;
            }
        }
    }

    xc_free_functional(xcfun);
    return 0;
}

void XCint::distribute_matrix2(const int block_length,
                               const int num_variables,
                               const int num_perturbations,
                               const int mat_dim,
                               const double prefactors[],
                               const int w_off,
                               const bool n_is_used[],
                               const double n[],
                               double u[],
                               double vxc[],
                               double &exc,
                               const std::vector<int> coor,
                               const double grid_w[])
//   const double grid_w[]) const
{
    int off;

    dens_offset = fun.set_order(num_perturbations + 1, xcfun);

    double *xcin = new double[num_variables * dens_offset * block_length];
    std::fill(&xcin[0], &xcin[num_variables * dens_offset * block_length], 0.0);

    double *xcout = new double[dens_offset * block_length];

    // has to be AO_BLOCK_LENGTH otherwise u can be too short
    std::fill(&u[0], &u[AO_BLOCK_LENGTH * num_variables], 0.0);

    for (int k = 0; k < MAX_NUM_DENSITIES; k++)
    {
        if (n_is_used[k])
        {
            for (int ivar = 0; ivar < num_variables; ivar++)
            {
                for (int ib = 0; ib < block_length; ib++)
                {
                    xcin[ib * num_variables * dens_offset + ivar * dens_offset +
                         k] = n[k * block_length * num_variables +
                                ivar * block_length + ib];
                }
            }
        }
    }

    for (int ivar = 0; ivar < num_variables; ivar++)
    {
        for (int jvar = 0; jvar < num_variables; jvar++)
        {
            off = jvar * dens_offset + (int)pow(2, num_perturbations);
            if (ivar == jvar)
            {
                for (int ib = 0; ib < block_length; ib++)
                {
                    xcin[off + ib * num_variables * dens_offset] = 1.0;
                }
            }
            else
            {
                for (int ib = 0; ib < block_length; ib++)
                {
                    xcin[off + ib * num_variables * dens_offset] = 0.0;
                }
            }
        }

        off = ivar * block_length;
        std::fill(&xcout[0], &xcout[dens_offset * block_length], 0.0);
        for (int ib = 0; ib < block_length; ib++)
        {
            if (n[ib] > 1.0e-14 and std::abs(grid_w[w_off + ib]) > 1.0e-30)
            {
                xc_eval(xcfun,
                        &xcin[ib * num_variables * dens_offset],
                        &xcout[ib * dens_offset]);
                u[off + ib] +=
                    xcout[(ib + 1) * dens_offset - 1] * grid_w[w_off + ib];
            }
        }
    }

    bool distribute_gradient;
    bool distribute_tau;

    if (fun.is_tau_mgga)
    {
        distribute_gradient = true;
        distribute_tau = true;
    }
    else if (fun.is_gga)
    {
        distribute_gradient = true;
        distribute_tau = false;
    }
    else
    {
        distribute_gradient = false;
        distribute_tau = false;
    }

    if (coor.size() == 0)
    {
        // FIXME can be moved one layer up since we need it also for density
        int max_ao_geo_order = 5; // FIXME hardcoded
        int buffer_len = balboa_get_buffer_len(
            balboa_context, max_ao_geo_order, AO_BLOCK_LENGTH);
        double *ao_compressed = new double[buffer_len];
        int *ao_compressed_index = new int[buffer_len];
        int ao_compressed_num;
        int num_aos = balboa_get_num_aos(balboa_context);
        int *ao_centers = new int[num_aos];
        for (int i = 0; i < num_aos; i++)
        {
            ao_centers[i] = balboa_get_ao_center(balboa_context, i);
        }
        int slice_offsets[4];
        compute_slice_offsets(std::vector<int>(), slice_offsets);
        compress(distribute_gradient,
                 block_length,
                 ao_compressed_num,
                 ao_compressed_index,
                 ao_compressed,
                 num_aos,
                 ao,
                 ao_centers,
                 std::vector<int>(),
                 slice_offsets);
        distribute_matrix(mat_dim,
                          block_length,
                          distribute_gradient,
                          distribute_tau,
                          prefactors,
                          u,
                          vxc,
                          ao_compressed_num,
                          ao_compressed_index,
                          ao_compressed,
                          ao_compressed_num,
                          ao_compressed_index,
                          ao_compressed);
        delete[] ao_compressed;
        delete[] ao_compressed_index;
        delete[] ao_centers;
    }
    else
    {
        auto get_geo_offset = [&](int i, int j, int k) {
            return balboa_get_geo_offset(balboa_context, i, j, k);
        };
        int num_aos = balboa_get_num_aos(balboa_context);
        int max_ao_geo_order = 5; // FIXME hardcoded
        int buffer_len = balboa_get_buffer_len(
            balboa_context, max_ao_geo_order, AO_BLOCK_LENGTH);
        int *ao_centers = new int[num_aos];
        for (int i = 0; i < num_aos; i++)
        {
            ao_centers[i] = balboa_get_ao_center(balboa_context, i);
        }
        get_mat_geo_derv(mat_dim,
                         num_aos,
                         block_length,
                         buffer_len,
                         ao,
                         ao_centers,
                         distribute_gradient,
                         distribute_tau,
                         coor,
                         get_geo_offset,
                         u,
                         vxc);
        delete[] ao_centers;
    }

    for (int ib = 0; ib < block_length; ib++)
    {
        exc += xcout[ib * dens_offset] * grid_w[w_off + ib];
    }

    delete[] xcin;
    delete[] xcout;
}

void XCint::compute_slice_offsets(const std::vector<int> &coor, int off[])
{
    int kp[3] = {0, 0, 0};
    for (size_t j = 0; j < coor.size(); j++)
    {
        kp[(coor[j] - 1) % 3]++;
    }
    off[0] = balboa_get_geo_offset(balboa_context, kp[0], kp[1], kp[2]);
    off[1] = balboa_get_geo_offset(balboa_context, kp[0] + 1, kp[1], kp[2]);
    off[2] = balboa_get_geo_offset(balboa_context, kp[0], kp[1] + 1, kp[2]);
    off[3] = balboa_get_geo_offset(balboa_context, kp[0], kp[1], kp[2] + 1);
}
