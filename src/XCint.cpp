#include "xcint.h"
#include "xcfun.h"

#include <math.h>
#include <time.h>
#include <assert.h>

#include <cstdlib>
#include <fstream>
#include <algorithm>

#include "rolex.h"
#include "XCint.h"
#include "AOBatch.h"
#include "MemAllocator.h"

#include "xcint_parameters.h"
#include "parameters.h"

#ifdef ENABLE_OMP
#include "omp.h"
#endif

#define AS_TYPE(Type, Obj) reinterpret_cast<Type *>(Obj)
#define AS_CTYPE(Type, Obj) reinterpret_cast<const Type *>(Obj)


// FIXME ugly global
xc_functional xcfun;
int dens_offset;


xcint_context_t *xcint_new()
{
    return AS_TYPE(xcint_context_t, new XCint());
}
XCint::XCint()
{
    nullify();
}


void xcint_free(xcint_context_t *context)
{
    if (!context) return;
    delete AS_TYPE(XCint, context);
}
XCint::~XCint()
{
    nullify();
}


XCINT_API int xcint_set_functional(xcint_context_t *context,
                                   const char      *line)
{
    return AS_TYPE(XCint, context)->set_functional(line);
}


XCINT_API int xcint_set_basis(xcint_context_t     *context,
                              const xcint_basis_t basis_type,
                              const int           num_centers,
                              const double        center_coordinates[],
                              const int           center_elements[],
                              const int           num_shells,
                              const int           shell_centers[],
                              const int           shell_l_quantum_numbers[],
                              const int           shell_num_primitives[],
                              const double        primitive_exponents[],
                              const double        contraction_coefficients[])
{
    return AS_TYPE(XCint, context)->set_basis(basis_type,
                                              num_centers,
                                              center_coordinates,
                                              center_elements,
                                              num_shells,
                                              shell_centers,
                                              shell_l_quantum_numbers,
                                              shell_num_primitives,
                                              primitive_exponents,
                                              contraction_coefficients);
}


void XCint::nullify()
{
    reset_time();
}


int XCint::set_functional(const char *line)
{
    fun.set_functional(line);
    return 0;
}


int XCint::set_basis(const int    basis_type,
                     const int    num_centers,
                     const double center_coordinates[],
                     const int    center_elements[],
                     const int    num_shells,
                     const int    shell_centers[],
                     const int    shell_l_quantum_numbers[],
                     const int    shell_num_primitives[],
                     const double primitive_exponents[],
                     const double contraction_coefficients[])
{
    basis.init(basis_type,
               num_centers,
               center_coordinates,
               center_elements,
               num_shells,
               shell_centers,
               shell_l_quantum_numbers,
               shell_num_primitives,
               primitive_exponents,
               contraction_coefficients);

    return 0;
}


void XCint::integrate_batch(const double dmat[],
                            const int    get_exc,
                                  double &exc,
                            const int    get_vxc,
                                  double vxc[],
                                  double &num_electrons,
                            const int    geo_coor[],
                            const bool   use_dmat[],
                            const int    ipoint,
                            const int    geo_derv_order,
                            const int    max_ao_order_g,
                            const int    block_length,
                            const int    num_variables,
                            const int    num_perturbations,
                            const int    num_fields,
                            const int    mat_dim,
                            const bool   get_gradient,
                            const bool   get_tau,
                            const int    dmat_index[],
                            const double grid[]) const
{
    AOBatch batch;

    size_t block_size = AO_BLOCK_LENGTH*num_variables*MAX_NUM_DENSITIES*sizeof(double);
    double *n = (double*) MemAllocator::allocate(block_size);
    block_size = AO_BLOCK_LENGTH*num_variables*sizeof(double);
    double *u = (double*) MemAllocator::allocate(block_size);

    std::vector<int> coor;
    double prefactors[5] = {1.0, 2.0, 2.0, 2.0, 0.5};

    bool n_is_used[MAX_NUM_DENSITIES];
    std::fill(&n_is_used[0], &n_is_used[MAX_NUM_DENSITIES], false);

    rolex::start_partial();

    batch.get_ao(basis,
                 get_gradient,
                 max_ao_order_g,
                 &grid[ipoint*4]);

//  time_ao += rolex::stop_partial();

    rolex::start_partial();

    if (!n_is_used[0])
    {
        std::fill(&n[0], &n[block_length*num_variables], 0.0);
        n_is_used[0] = true;
    }

    batch.get_density_undiff(mat_dim,
                             get_gradient,
                             get_tau,
                             prefactors,
                             n,
                             dmat,
                             true,
                             true);

    for (int ib = 0; ib < block_length; ib++)
    {
        num_electrons += grid[(ipoint + ib)*4 + 3]*n[ib];
    }

//  time_densities += rolex::stop_partial();

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
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                batch.get_density_undiff(mat_dim,
                                         get_gradient,
                                         get_tau,
                                         prefactors,
                                         &n[k*block_length*num_variables],
                                         &dmat[dmat_index[k]],
                                         false,
                                         false); // FIXME can be true based on dmat, saving possible
            }
        }

        rolex::start_partial();

#include "ave_contributions.h"

//      time_densities += rolex::stop_partial();

        dens_offset = fun.set_order(num_perturbations, xcfun);

        size_t block_size = num_variables*dens_offset*block_length*sizeof(double);
        xcin = (double*) MemAllocator::allocate(block_size);
        std::fill(&xcin[0], &xcin[num_variables*dens_offset*block_length], 0.0);

        block_size = dens_offset*block_length*sizeof(double);
        xcout = (double*) MemAllocator::allocate(block_size);

        for (int k = 0; k < MAX_NUM_DENSITIES; k++)
        {
            if (n_is_used[k])
            {
                for (int ivar = 0; ivar < num_variables; ivar++)
                {
                    for (int ib = 0; ib < block_length; ib++)
                    {
                        xcin[ib*num_variables*dens_offset + ivar*dens_offset + k] = n[k*block_length*num_variables + ivar*block_length + ib];
                    }
                }
            }
        }

        rolex::start_partial();

        double sum = 0.0;
        for (int ib = 0; ib < block_length; ib++)
        {
            if (n[ib] > 1.0e-14 and fabs(grid[(ipoint + ib)*4 + 3]) > 1.0e-30)
            {
                xc_eval(xcfun, &xcin[ib*num_variables*dens_offset], &xcout[ib*dens_offset]);
                sum += xcout[ib*dens_offset + dens_offset - 1]*grid[(ipoint + ib)*4 + 3];
            }
        }
        exc += sum;

//      time_fun_derv += rolex::stop_partial();

        MemAllocator::deallocate(xcin);
        MemAllocator::deallocate(xcout);
    }


    // matrix contribution
    if (get_vxc)
    {
        int k;

        bool contribution_is_implemented = false;

        if (geo_derv_order == 0) // no geo dervs
        {
            contribution_is_implemented = true;

            rolex::start_partial();

            for (int ifield = 0; ifield < num_fields; ifield++)
            {
                k = (int)pow(2, ifield); // FIXME double check this
                if (!n_is_used[k])
                {
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                batch.get_density_undiff(mat_dim,
                                         get_gradient,
                                         get_tau,
                                         prefactors,
                                         &n[k*block_length*num_variables],
                                         &dmat[(ifield+1)*mat_dim*mat_dim],
                                         false,
                                         false); // FIXME can be true depending on perturbation (savings possible)
            }

//          time_densities += rolex::stop_partial();

            distribute_matrix(block_length,
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
                              batch,
                              grid);
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
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                batch.get_dens_geo_derv(basis,
                                        mat_dim,
                                        get_gradient,
                                        get_tau,
                                        coor,
                                        &n[k*block_length*num_variables],
                                        &dmat[0]);
                coor.clear();
                distribute_matrix(block_length,
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
                                  batch,
                                  grid);
                n_is_used[1] = false;

                // M_i  d_n
                coor.push_back(geo_coor[0]);
                distribute_matrix(block_length,
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
                                  batch,
                                  grid);
                coor.clear();
            }

            if (geo_derv_order == 2 && num_fields == 0)
            {
                contribution_is_implemented = true;

                // M_ij d_n
                coor.push_back(geo_coor[0]);
                coor.push_back(geo_coor[1]);
                distribute_matrix(block_length,
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
                                  batch,
                                  grid);
                coor.clear();

                // FIXME add shortcut if i == j
                // M_i  d_nn  n_j
                k = 1;
                if (!n_is_used[k])
                {
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[1]);
                batch.get_dens_geo_derv(basis,
                                        mat_dim,
                                        get_gradient,
                                        get_tau,
                                        coor,
                                        &n[k*block_length*num_variables],
                                        &dmat[0]);
                coor.clear();
                coor.push_back(geo_coor[0]);
                distribute_matrix(block_length,
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
                                  batch,
                                  grid);
                coor.clear();
                n_is_used[1] = false;

                // M_j  d_nn  n_i
                k = 1;
                if (!n_is_used[k])
                {
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                batch.get_dens_geo_derv(basis,
                                        mat_dim,
                                        get_gradient,
                                        get_tau,
                                        coor,
                                        &n[k*block_length*num_variables],
                                        &dmat[0]);
                coor.clear();
                coor.push_back(geo_coor[1]);
                distribute_matrix(block_length,
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
                                  batch,
                                  grid);
                coor.clear();
                n_is_used[1] = false;

                // M    d_nnn n_i n_j
                // M    d_nn  n_ij
                // FIXME stupidly recalculating n_i and n_j
                k = 1;
                if (!n_is_used[k])
                {
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                batch.get_dens_geo_derv(basis,
                                        mat_dim,
                                        get_gradient,
                                        get_tau,
                                        coor,
                                        &n[k*block_length*num_variables],
                                        &dmat[0]);
                coor.clear();
                k = 2;
                if (!n_is_used[k])
                {
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[1]);
                batch.get_dens_geo_derv(basis,
                                        mat_dim,
                                        get_gradient,
                                        get_tau,
                                        coor,
                                        &n[k*block_length*num_variables],
                                        &dmat[0]);
                coor.clear();
                k = 3;
                if (!n_is_used[k])
                {
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                coor.push_back(geo_coor[1]);
                batch.get_dens_geo_derv(basis,
                                        mat_dim,
                                        get_gradient,
                                        get_tau,
                                        coor,
                                        &n[k*block_length*num_variables],
                                        &dmat[0]);
                coor.clear();
                distribute_matrix(block_length,
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
                                  batch,
                                  grid);
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
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                batch.get_density_undiff(mat_dim,
                                         get_gradient,
                                         get_tau,
                                         prefactors,
                                         &n[k*block_length*num_variables],
                                         &dmat[1*mat_dim*mat_dim],
                                         false,
                                         false); // FIXME can be true depending on perturbation (savings possible)
                coor.push_back(geo_coor[0]);
                distribute_matrix(block_length,
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
                                  batch,
                                  grid);
                coor.clear();

                // M    d_nnn n_i n_a
                // M    d_nn  n_ia
                k = 2;
                if (!n_is_used[k])
                {
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                batch.get_dens_geo_derv(basis,
                                        mat_dim,
                                        get_gradient,
                                        get_tau,
                                        coor,
                                        &n[k*block_length*num_variables],
                                        &dmat[0]);
                coor.clear();
                k = 3;
                if (!n_is_used[k])
                {
                    std::fill(&n[k*block_length*num_variables], &n[(k+1)*block_length*num_variables], 0.0);
                    n_is_used[k] = true;
                }
                coor.push_back(geo_coor[0]);
                batch.get_dens_geo_derv(basis,
                                        mat_dim,
                                        get_gradient,
                                        get_tau,
                                        coor,
                                        &n[k*block_length*num_variables],
                                        &dmat[1*mat_dim*mat_dim]);
                coor.clear();
                distribute_matrix(block_length,
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
                                  batch,
                                  grid);
                n_is_used[1] = false;
                n_is_used[2] = false;
                n_is_used[3] = false;
            }
        }

        if (!contribution_is_implemented)
        {
            fprintf(stderr, "ERROR: XCint matrix contribution for geo_derv_order=%i and num_fields=%i not implemented\n", geo_derv_order, num_fields);
            exit(-1);
        }
    }

    MemAllocator::deallocate(n);
    MemAllocator::deallocate(u);
}


XCINT_API int xcint_integrate(const xcint_context_t      *context,
                              const xcint_mode_t         mode,
                              const int                  num_points,
                              const double               grid[],
                              const int                  num_perturbations,
                              const xcint_perturbation_t perturbations[],
                              const int                  components[],
                              const int                  num_dmat,
                              const int                  dmat_to_perturbations[],
                              const int                  dmat_to_components[],
                              const double               dmat[],
                              const int                  get_exc,
                                    double               *exc,
                              const int                  get_vxc,
                                    double               vxc[],
                                    double               *num_electrons)
{
    return AS_CTYPE(XCint, context)->integrate(mode,
                                               num_points,
                                               grid,
                                               num_perturbations,
                                               perturbations,
                                               components,
                                               num_dmat,
                                               dmat_to_perturbations,
                                               dmat_to_components,
                                               dmat,
                                               get_exc,
                                               exc,
                                               get_vxc,
                                               vxc,
                                               num_electrons);
}
int XCint::integrate(const xcint_mode_t         mode,
                     const int                  num_points,
                     const double               grid[],
                     const int                  num_perturbations,
                     const xcint_perturbation_t perturbations[],
                     const int                  components[],
                     const int                  num_dmat,
                     const int                  dmat_to_perturbations[],
                     const int                  dmat_to_components[],
                     const double               dmat[],
                     const int                  get_exc,
                           double               *exc,
                     const int                  get_vxc,
                           double               vxc[],
                           double               *num_electrons) const
{
    xcfun = xc_new_functional();
    for (int i = 0; i < fun.keys.size(); i++)
    {
        int ierr = xc_set(xcfun, fun.keys[i].c_str(), fun.weights[i]);
        if (ierr != 0)
        {
            fprintf(stderr, "ERROR in fun init: \"%s\" not recognized, quitting.\n", fun.keys[i].c_str());
            exit(-1);
        }
    }

    int num_variables;
    int max_ao_order_g;
    size_t block_size;

    double a;
    double prefactors[5] = {1.0, 2.0, 2.0, 2.0, 0.5};

    double *n = NULL;
    double *u = NULL;

    bool contribution_is_implemented;

    std::vector<int> coor;

    rolex::start_global();

    int rank = 0;
    int num_proc = 1;

    int mat_dim;
    if (rank == 0) mat_dim = basis.get_num_ao();

    int geo_derv_order = 0;
    int num_fields = 0;
    for (int i = 0; i < num_perturbations; i++)
    {
        if (perturbations[i] == XCINT_PERT_GEO) geo_derv_order++;
        if (perturbations[i] == XCINT_PERT_EL)  num_fields++;
    }

    int *geo_coor = NULL;
    if (geo_derv_order > 0)
    {
        block_size = geo_derv_order*sizeof(int);
        geo_coor = (int*) MemAllocator::allocate(block_size);
        for (int i = 0; i < num_perturbations; i++)
        {
            if (perturbations[i] == XCINT_PERT_GEO) geo_coor[i] = components[2*i]; // FIXME
        }
    }

    assert(mode == XCINT_MODE_RKS);

    *exc = 0.0;
    *num_electrons = 0.0;

    if (get_vxc) std::fill(&vxc[0], &vxc[mat_dim*mat_dim], 0.0);

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

//  FIXME
//  this causes the geo_off array to become
//  too small and leads to out of bounds access
//  basis.set_geo_off(max_ao_order_g);

    rolex::start_partial();

    bool *use_dmat = NULL;
    int  *dmat_index = NULL;

    block_size = MAX_NUM_DENSITIES*sizeof(int);

    dmat_index = (int*) MemAllocator::allocate(block_size);
    std::fill(&dmat_index[0], &dmat_index[MAX_NUM_DENSITIES], 0);

    use_dmat = (bool*) MemAllocator::allocate(block_size);
    std::fill(&use_dmat[0], &use_dmat[MAX_NUM_DENSITIES], false);

    assert(num_dmat <= MAX_NUM_DENSITIES);

    for (int k = 0; k < num_dmat; k++)
    {
        use_dmat[dmat_to_perturbations[k]] = true;
        dmat_index[dmat_to_perturbations[k]] = k*mat_dim*mat_dim;
    }

    assert(num_perturbations < 7);

    int block_length;

#ifdef ENABLE_OMP
    int num_threads = 0;

    #pragma omp parallel
    {
        if (omp_get_thread_num() == 0) num_threads = omp_get_num_threads();
    }

    double *exc_buffer = NULL;
    if (get_exc) exc_buffer = (double*) MemAllocator::allocate(num_threads*sizeof(double));

    double *num_electrons_buffer = (double*) MemAllocator::allocate(num_threads*sizeof(double));

    double *vxc_buffer = NULL;
    if (get_vxc)
    {
        vxc_buffer = (double*) MemAllocator::allocate(num_threads*mat_dim*mat_dim*sizeof(double));
        std::fill(&vxc_buffer[0], &vxc_buffer[num_threads*mat_dim*mat_dim], 0.0);
    }

    #pragma omp parallel
    {
        int ithread = omp_get_thread_num();

        double exc_local = 0.0;
        double num_electrons_local = 0.0;

        double *vxc_local = NULL;
        if (get_vxc) vxc_local = &vxc_buffer[ithread*mat_dim*mat_dim];

        #pragma omp for schedule(dynamic)
#else
        double exc_local = *exc;
        double num_electrons_local = *num_electrons;
        double *vxc_local = NULL;
        if (get_vxc) vxc_local = &vxc[0];
#endif
        for (int ibatch = 0; ibatch < num_points/AO_BLOCK_LENGTH; ibatch++)
        {
            int ipoint = ibatch*AO_BLOCK_LENGTH;

            block_length = AO_BLOCK_LENGTH;

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
                            grid);
        }

#ifdef ENABLE_OMP
        if (get_exc) exc_buffer[ithread] = exc_local;
        num_electrons_buffer[ithread] = num_electrons_local;
    }

    for (size_t ithread = 0; ithread < num_threads; ithread++)
    {
        num_electrons += num_electrons_buffer[ithread];
        if (get_exc) exc += exc_buffer[ithread];
        if (get_vxc)
        {
            // FIXME consider using blas daxpy for this
            for (int i = 0; i < mat_dim*mat_dim; i++)
            {
                vxc[i] += vxc_buffer[ithread*mat_dim*mat_dim + i];
            }
        }
    }

    MemAllocator::deallocate(num_electrons_buffer);
    MemAllocator::deallocate(exc_buffer);
    MemAllocator::deallocate(vxc_buffer);
#else
    *exc = exc_local;
    *num_electrons = num_electrons_local;
#endif

    MemAllocator::deallocate(use_dmat);
    MemAllocator::deallocate(dmat_index);
    MemAllocator::deallocate(geo_coor);

    if (get_vxc)
    {
        // symmetrize result matrix
        if (rank == 0)
        {
            for (int k = 0; k < mat_dim; k++)
            {
                for (int l = 0; l < k; l++)
                {
                    a = vxc[k*mat_dim + l] + vxc[l*mat_dim + k];
                    vxc[k*mat_dim + l] = 0.5*a;
                    vxc[l*mat_dim + k] = 0.5*a;
                }
            }
        }
    }

//  time_total += rolex::stop_global();
    xc_free_functional(xcfun);
    return 0;
}


void XCint::distribute_matrix(const int              block_length,
                              const int              num_variables,
                              const int              num_perturbations,
                              const int              mat_dim,
                              const double           prefactors[],
                              const int              w_off,
                              const bool             n_is_used[],
                              const double           n[],
                                    double           u[],
                                    double           vxc[],
                                    double           &exc,
                              const std::vector<int> coor,
                                    AOBatch     &batch,
                              const double           grid[]) const
{
    rolex::start_partial();

    int off;

    double *xcin = NULL;
    double *xcout = NULL;
    size_t block_size;

    dens_offset = fun.set_order(num_perturbations + 1, xcfun);

    block_size = num_variables*dens_offset*block_length*sizeof(double);
    xcin = (double*) MemAllocator::allocate(block_size);
    std::fill(&xcin[0], &xcin[num_variables*dens_offset*block_length], 0.0);

    block_size = dens_offset*block_length*sizeof(double);
    xcout = (double*) MemAllocator::allocate(block_size);

    std::fill(&u[0], &u[block_length*num_variables], 0.0);

    for (int k = 0; k < MAX_NUM_DENSITIES; k++)
    {
        if (n_is_used[k])
        {
            for (int ivar = 0; ivar < num_variables; ivar++)
            {
                for (int ib = 0; ib < block_length; ib++)
                {
                    xcin[ib*num_variables*dens_offset + ivar*dens_offset + k] = n[k*block_length*num_variables + ivar*block_length + ib];
                }
            }
        }
    }

    for (int ivar = 0; ivar < num_variables; ivar++)
    {
        for (int jvar = 0; jvar < num_variables; jvar++)
        {
            off = jvar*dens_offset + (int)pow(2, num_perturbations);
            if (ivar == jvar)
            {
                for (int ib = 0; ib < block_length; ib++)
                {
                    xcin[off + ib*num_variables*dens_offset] = 1.0;
                }
            }
            else
            {
                for (int ib = 0; ib < block_length; ib++)
                {
                    xcin[off + ib*num_variables*dens_offset] = 0.0;
                }
            }
        }

        off = ivar*block_length;
        std::fill(&xcout[0], &xcout[dens_offset*block_length], 0.0);
        for (int ib = 0; ib < block_length; ib++)
        {
            if (n[ib] > 1.0e-14 and fabs(grid[(w_off + ib)*4 + 3]) > 1.0e-30)
            {
                xc_eval(xcfun, &xcin[ib*num_variables*dens_offset], &xcout[ib*dens_offset]);
                u[off + ib] += xcout[(ib+1)*dens_offset - 1]*grid[(w_off + ib)*4 + 3];
            }
        }
    }

//  time_fun_derv += rolex::stop_partial();

    rolex::start_partial();

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
        batch.distribute_matrix_undiff(mat_dim,
                                       distribute_gradient,
                                       distribute_tau,
                                       prefactors,
                                       u,
                                       vxc);
    }
    else
    {
        batch.get_mat_geo_derv(basis,
                               mat_dim,
                               distribute_gradient,
                               distribute_tau,
                               coor,
                               u,
                               vxc);
    }

    for (int ib = 0; ib < block_length; ib++)
    {
        exc += xcout[ib*dens_offset]*grid[(w_off + ib)*4 + 3];
    }

    MemAllocator::deallocate(xcin);
    MemAllocator::deallocate(xcout);

//  time_matrix_distribution += rolex::stop_partial();
}


void XCint::reset_time()
{
//  time_total = 0.0;
//  time_ao = 0.0;
//  time_fun_derv = 0.0;
//  time_densities = 0.0;
//  time_matrix_distribution = 0.0;
}
