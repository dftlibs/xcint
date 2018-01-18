#include <math.h>
#include <stdio.h>

#include <algorithm>
#include <assert.h>
#include <iostream>

#include "AOBatch.h"
#include "blas_interface.h"
#include "density_parameters.h"
#include "compress.h"

AOBatch::AOBatch()
{
    ao_length = -1;
    ao = NULL;
    assert(AO_BLOCK_LENGTH % AO_CHUNK_LENGTH == 0);
    balboa_context = balboa_new_context();
}

AOBatch::~AOBatch()
{
    delete[] ao;
    ao_length = -1;
    ao = NULL;
    balboa_free_context(balboa_context);
}

void AOBatch::get_ao(const bool use_gradient,
                     const int max_ao_geo_order,
                     const int block_length,
                     const double grid_x_bohr[],
                     const double grid_y_bohr[],
                     const double grid_z_bohr[])
{
    assert(max_ao_geo_order <= MAX_GEO_DIFF_ORDER);

    int buffer_len = balboa_get_buffer_len(
        balboa_context, max_ao_geo_order, AO_BLOCK_LENGTH);
    //  FIXME should be:
    //  int buffer_len = balboa_get_buffer_len(balboa_context, max_ao_geo_order,
    //  block_length);

    if (buffer_len != ao_length) // FIXME
    {
        ao_length = buffer_len;
        delete[] ao;
        ao = new double[buffer_len];
    }

    std::fill(&ao[0], &ao[buffer_len], 0.0);

    // FIXME is this really needed?
    double *buffer = new double[buffer_len];
    std::fill(&buffer[0], &buffer[buffer_len], 0.0);

    int ierr;
    ierr = balboa_get_ao(balboa_context,
                         max_ao_geo_order,
                         block_length,
                         grid_x_bohr,
                         grid_y_bohr,
                         grid_z_bohr,
                         buffer);

    std::copy(&buffer[0], &buffer[buffer_len], &ao[0]);

    delete[] buffer;
}

void AOBatch::distribute_matrix_undiff(const int mat_dim,
                                       const bool use_gradient,
                                       const bool use_tau,
                                       const double prefactors[],
                                       const double u[],
                                       double fmat[])
{
    // FIXME can be moved one layer up since we need it also for density
    int max_ao_geo_order = 5; // FIXME hardcoded
    int buffer_len = balboa_get_buffer_len(balboa_context, max_ao_geo_order, AO_BLOCK_LENGTH);
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
    compress(use_gradient,
             ao_compressed_num,
             ao_compressed_index,
             ao_compressed,
             num_aos,
             ao,
             ao_centers,
             std::vector<int>(),
             slice_offsets);
    distribute_matrix(mat_dim,
                      use_gradient,
                      use_tau,
                      prefactors,
                      u,
                      fmat,
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

void AOBatch::distribute_matrix(const int mat_dim,
                                const bool use_gradient,
                                const bool use_tau,
                                const double prefactors[],
                                const double u[],
                                double fmat[],
                                const int k_aoc_num,
                                const int k_aoc_index[],
                                const double k_aoc[],
                                const int l_aoc_num,
                                const int l_aoc_index[],
                                const double l_aoc[])
{
    // here we compute       F(k, l) += AO_k(k, b) u(b) AO_l(l, b)
    // in two steps
    // step 1:               W(k, b)  = AO_k(k, b) u(b)
    // step 2:               F(k, l) += W(k, b) AO_l(l, b)^T

    if (k_aoc_num == 0)
        return;
    if (l_aoc_num == 0)
        return;

    double *W = new double[k_aoc_num * AO_BLOCK_LENGTH];

    std::fill(&W[0], &W[AO_BLOCK_LENGTH * k_aoc_num], 0.0);

    int kc, lc, iboff;

    int num_slices;
    (use_gradient) ? (num_slices = 4) : (num_slices = 1);

    for (int islice = 0; islice < num_slices; islice++)
    {
        if (fabs(prefactors[islice]) > 0.0)
        {
            for (int k = 0; k < k_aoc_num; k++)
            {
                iboff = k * AO_BLOCK_LENGTH;
                for (int ib = 0; ib < AO_BLOCK_LENGTH; ib++)
                {
                    W[iboff + ib] +=
                        prefactors[islice] * u[islice * AO_BLOCK_LENGTH + ib] *
                        k_aoc[islice * AO_BLOCK_LENGTH * mat_dim + iboff + ib];
                }
            }
        }
    }

    double *F = new double[k_aoc_num * l_aoc_num];

    // we compute F(k, l) += W(k, b) AO_l(l, b)^T
    // we transpose W instead of AO_l because we call fortran blas
    char ta = 't';
    char tb = 'n';
    int im = k_aoc_num;
    int in = l_aoc_num;
    int ik = AO_BLOCK_LENGTH;
    int lda = ik;
    int ldb = ik;
    int ldc = im;
    double alpha = 1.0;
    double beta = 0.0;
    wrap_dgemm(
        &ta, &tb, &im, &in, &ik, &alpha, W, &lda, l_aoc, &ldb, &beta, F, &ldc);

    if (use_tau)
    {
        if (fabs(prefactors[4]) > 0.0)
        {
            for (int ixyz = 0; ixyz < 3; ixyz++)
            {
                for (int k = 0; k < k_aoc_num; k++)
                {
                    iboff = k * AO_BLOCK_LENGTH;
                    for (int ib = 0; ib < AO_BLOCK_LENGTH; ib++)
                    {
                        W[iboff + ib] =
                            u[4 * AO_BLOCK_LENGTH + ib] *
                            k_aoc[(ixyz + 1) * AO_BLOCK_LENGTH * mat_dim +
                                  iboff + ib];
                    }
                }

                alpha = prefactors[4];
                beta = 1.0;
                wrap_dgemm(&ta,
                           &tb,
                           &im,
                           &in,
                           &ik,
                           &alpha,
                           W,
                           &lda,
                           &l_aoc[(ixyz + 1) * AO_BLOCK_LENGTH * mat_dim],
                           &ldb,
                           &beta,
                           F,
                           &ldc);
            }
        }
    }

    delete[] W;
    W = NULL;

    // FIXME easier route possible if k and l match
    for (int k = 0; k < k_aoc_num; k++)
    {
        kc = k_aoc_index[k];
        for (int l = 0; l < l_aoc_num; l++)
        {
            lc = l_aoc_index[l];
            fmat[kc * mat_dim + lc] += F[k * l_aoc_num + l];
        }
    }

    delete[] F;
    F = NULL;
}

void AOBatch::get_density_undiff(const int mat_dim,
                                 const bool use_gradient,
                                 const bool use_tau,
                                 const double prefactors[],
                                 double density[],
                                 const double dmat[],
                                 const bool dmat_is_symmetric,
                                 const bool kl_match)
{
    // FIXME can be moved one layer up since we need it also for matrix distribution
    int max_ao_geo_order = 5; // FIXME hardcoded
    int buffer_len = balboa_get_buffer_len(balboa_context, max_ao_geo_order, AO_BLOCK_LENGTH);
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
    compress(use_gradient,
             ao_compressed_num,
             ao_compressed_index,
             ao_compressed,
             num_aos,
             ao,
             ao_centers,
             std::vector<int>(),
             slice_offsets);
    get_density(mat_dim,
                use_gradient,
                use_tau,
                prefactors,
                density,
                dmat,
                dmat_is_symmetric,
                kl_match,
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

void AOBatch::get_density(const int mat_dim,
                          const bool use_gradient,
                          const bool use_tau,
                          const double prefactors[],
                          double density[],
                          const double dmat[],
                          const bool dmat_is_symmetric,
                          const bool kl_match,
                          const int k_aoc_num,
                          const int k_aoc_index[],
                          const double k_aoc[],
                          const int l_aoc_num,
                          const int l_aoc_index[],
                          const double l_aoc[])
{
    // here we compute       n(b)    = AO_k(k, b) D(k, l) AO_l(l, b)
    // in two steps
    // step 1:               X(k, b) = D(k, l) AO_l(l, b)
    // step 2:               n(b)    = AO_k(k, b) X(k, b)

    if (k_aoc_num == 0)
        return;
    if (l_aoc_num == 0)
        return;

    int kc, lc, iboff;

    double *D = new double[k_aoc_num * l_aoc_num];
    double *X = new double[k_aoc_num * AO_BLOCK_LENGTH];

    // compress dmat
    if (kl_match)
    {
        if (dmat_is_symmetric)
        {
            for (int k = 0; k < k_aoc_num; k++)
            {
                kc = k_aoc_index[k];
                for (int l = 0; l <= k; l++)
                {
                    lc = l_aoc_index[l];
                    D[k * l_aoc_num + l] = 2.0 * dmat[kc * mat_dim + lc];
                }
            }
        }
        else
        {
            for (int k = 0; k < k_aoc_num; k++)
            {
                kc = k_aoc_index[k];
                for (int l = 0; l < l_aoc_num; l++)
                {
                    lc = l_aoc_index[l];
                    D[k * l_aoc_num + l] = 2.0 * dmat[kc * mat_dim + lc];
                }
            }
        }
    }
    else
    {
        if (dmat_is_symmetric)
        {
            for (int k = 0; k < k_aoc_num; k++)
            {
                kc = k_aoc_index[k];
                for (int l = 0; l < l_aoc_num; l++)
                {
                    lc = l_aoc_index[l];
                    D[k * l_aoc_num + l] = 2.0 * dmat[kc * mat_dim + lc];
                }
            }
        }
        else
        {
            for (int k = 0; k < k_aoc_num; k++)
            {
                kc = k_aoc_index[k];
                for (int l = 0; l < l_aoc_num; l++)
                {
                    lc = l_aoc_index[l];
                    D[k * l_aoc_num + l] =
                        dmat[kc * mat_dim + lc] + dmat[lc * mat_dim + kc];
                }
            }
        }
    }

    // form xmat
    if (dmat_is_symmetric)
    {
        char si = 'r';
        char up = 'u';
        int im = AO_BLOCK_LENGTH;
        int in = k_aoc_num;
        int lda = in;
        int ldb = im;
        int ldc = im;
        double alpha = 1.0;
        double beta = 0.0;
        wrap_dsymm(
            &si, &up, &im, &in, &alpha, D, &lda, l_aoc, &ldb, &beta, X, &ldc);
    }
    else
    {
        char ta = 'n';
        char tb = 'n';
        int im = AO_BLOCK_LENGTH;
        int in = k_aoc_num;
        int ik = l_aoc_num;
        int lda = im;
        int ldb = ik;
        int ldc = im;
        double alpha = 1.0;
        double beta = 0.0;
        wrap_dgemm(&ta,
                   &tb,
                   &im,
                   &in,
                   &ik,
                   &alpha,
                   l_aoc,
                   &lda,
                   D,
                   &ldb,
                   &beta,
                   X,
                   &ldc);
    }

    int num_slices;
    (use_gradient) ? (num_slices = 4) : (num_slices = 1);

    std::fill(&density[0], &density[num_slices*AO_BLOCK_LENGTH], 0.0);

    // assemble density and possibly gradient
    for (int islice = 0; islice < num_slices; islice++)
    {
        if (fabs(prefactors[islice]) > 0.0)
        {
            for (int k = 0; k < k_aoc_num; k++)
            {
                iboff = k * AO_BLOCK_LENGTH;
                for (int ib = 0; ib < AO_BLOCK_LENGTH; ib++)
                {
                    density[islice * AO_BLOCK_LENGTH + ib] +=
                        prefactors[islice] * X[iboff + ib] *
                        k_aoc[islice * AO_BLOCK_LENGTH * mat_dim + iboff + ib];
                }
            }
        }
    }

    if (use_tau)
    {
        if (fabs(prefactors[4]) > 0.0)
        {
            for (int ixyz = 0; ixyz < 3; ixyz++)
            {
                if (dmat_is_symmetric)
                {
                    char si = 'r';
                    char up = 'u';
                    int im = AO_BLOCK_LENGTH;
                    int in = k_aoc_num;
                    int lda = in;
                    int ldb = im;
                    int ldc = im;
                    double alpha = 1.0;
                    double beta = 0.0;
                    wrap_dsymm(&si,
                               &up,
                               &im,
                               &in,
                               &alpha,
                               D,
                               &lda,
                               &l_aoc[(ixyz + 1) * AO_BLOCK_LENGTH * mat_dim],
                               &ldb,
                               &beta,
                               X,
                               &ldc);
                }
                else
                {
                    char ta = 'n';
                    char tb = 'n';
                    int im = AO_BLOCK_LENGTH;
                    int in = k_aoc_num;
                    int ik = l_aoc_num;
                    int lda = im;
                    int ldb = ik;
                    int ldc = im;
                    double alpha = 1.0;
                    double beta = 0.0;
                    wrap_dgemm(&ta,
                               &tb,
                               &im,
                               &in,
                               &ik,
                               &alpha,
                               &l_aoc[(ixyz + 1) * AO_BLOCK_LENGTH * mat_dim],
                               &lda,
                               D,
                               &ldb,
                               &beta,
                               X,
                               &ldc);
                }

                for (int k = 0; k < k_aoc_num; k++)
                {
                    iboff = k * AO_BLOCK_LENGTH;
                    for (int ib = 0; ib < AO_BLOCK_LENGTH; ib++)
                    {
                        density[4 * AO_BLOCK_LENGTH + ib] +=
                            prefactors[4] * X[iboff + ib] *
                            k_aoc[(ixyz + 1) * AO_BLOCK_LENGTH * mat_dim +
                                  iboff + ib];
                    }
                }
            }
        }
    }

    delete[] D;
    D = NULL;
    delete[] X;
    X = NULL;
}

void AOBatch::get_dens_geo_derv(const int mat_dim,
                                const bool use_gradient,
                                const bool use_tau,
                                const std::vector<int> &coor,
                                double density[],
                                const double mat[])
{
    /*
    1st                        a,0
    2nd            ab,0                    a,b
                  /    \                  /    \
                 /      \                /      \
    3rd     abc,0        ab,c        ac,b        a,bc
           /    \       /   \       /   \       /   \
    4th abcd,0 abc,d abd,c ab,cd acd,b ac,bd ad,bc a,bcd
    */

    for (unsigned int i = 0; i < coor.size(); i++)
    {
        if (coor[i] == 0)
            return;
    }

    // there is a factor 2 because center-a
    // differentiation can be left or right
    double f = 2.0 * pow(-1.0, (int)coor.size());

    std::vector<int> k_coor;
    std::vector<int> l_coor;

    switch (coor.size())
    {
    // FIXME implement shortcuts
    case 1:
        k_coor.push_back(coor[0]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 2:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 3:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[1]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 4:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        k_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[1]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        break;
    default:
        std::cout << "ERROR: get_dens_geo_derv coor.size() too high\n";
        exit(1);
        break;
    }
}

void AOBatch::get_mat_geo_derv(const int mat_dim,
                               const bool use_gradient,
                               const bool use_tau,
                               const std::vector<int> &coor,
                               const double density[],
                               double mat[])
{
    /*
    1st                        a,0
    2nd            ab,0                    a,b
                  /    \                  /    \
                 /      \                /      \
    3rd     abc,0        ab,c        ac,b        a,bc
           /    \       /   \       /   \       /   \
    4th abcd,0 abc,d abd,c ab,cd acd,b ac,bd ad,bc a,bcd
    */

    for (unsigned int i = 0; i < coor.size(); i++)
    {
        if (coor[i] == 0)
            return;
    }

    // there is a factor 2 because center-a
    // differentiation can be left or right
    double f = 2.0 * pow(-1.0, (int)coor.size());

    std::vector<int> k_coor;
    std::vector<int> l_coor;

    switch (coor.size())
    {
    // FIXME implement shortcuts
    case 1:
        k_coor.push_back(coor[0]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 2:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 3:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[1]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 4:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        k_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[1]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(
            mat_dim, use_gradient, use_tau, f, k_coor, l_coor, density, mat);
        k_coor.clear();
        l_coor.clear();
        break;
    default:
        std::cout << "ERROR: get_dens_geo_derv coor.size() too high\n";
        exit(1);
        break;
    }
}

void AOBatch::diff_u_wrt_center_tuple(const int mat_dim,
                                      const bool use_gradient,
                                      const bool use_tau,
                                      const double f,
                                      const std::vector<int> &k_coor,
                                      const std::vector<int> &l_coor,
                                      double u[],
                                      const double M[])
{
    int max_ao_geo_order = 5; // FIXME hardcoded
    int buffer_len = balboa_get_buffer_len(balboa_context, max_ao_geo_order, AO_BLOCK_LENGTH);
    double *k_ao_compressed = new double[buffer_len];
    int *k_ao_compressed_index = new int[buffer_len];
    int k_ao_compressed_num;
    double *l_ao_compressed = new double[buffer_len];
    int *l_ao_compressed_index = new int[buffer_len];
    int l_ao_compressed_num;
    int num_aos = balboa_get_num_aos(balboa_context);
    int *ao_centers = new int[num_aos];
    for (int i = 0; i < num_aos; i++)
    {
        ao_centers[i] = balboa_get_ao_center(balboa_context, i);
    }
    int slice_offsets[4];
    compute_slice_offsets(k_coor, slice_offsets);
    compress(use_gradient,
             k_ao_compressed_num,
             k_ao_compressed_index,
             k_ao_compressed,
             num_aos,
             ao,
             ao_centers,
             k_coor,
             slice_offsets);
    compute_slice_offsets(l_coor, slice_offsets);
    compress(use_gradient,
             l_ao_compressed_num,
             l_ao_compressed_index,
             l_ao_compressed,
             num_aos,
             ao,
             ao_centers,
             l_coor,
             slice_offsets);

    double prefactors[5] = {f, f, f, f, 0.5 * f};

    // evaluate density dervs
    get_density(mat_dim,
                use_gradient,
                use_tau,
                prefactors,
                u,
                M,
                false,
                false,
                k_ao_compressed_num,
                k_ao_compressed_index,
                k_ao_compressed,
                l_ao_compressed_num,
                l_ao_compressed_index,
                l_ao_compressed);
    if (use_gradient)
    {
        prefactors[0] = 0.0;
        get_density(mat_dim,
                    use_gradient,
                    false,
                    prefactors,
                    u,
                    M,
                    false,
                    false,
                    l_ao_compressed_num,
                    l_ao_compressed_index,
                    l_ao_compressed,
                    k_ao_compressed_num,
                    k_ao_compressed_index,
                    k_ao_compressed);
    }
    delete[] k_ao_compressed;
    delete[] k_ao_compressed_index;
    delete[] l_ao_compressed;
    delete[] l_ao_compressed_index;
    delete[] ao_centers;
}

void AOBatch::diff_M_wrt_center_tuple(const int mat_dim,
                                      const bool use_gradient,
                                      const bool use_tau,
                                      const double f,
                                      const std::vector<int> &k_coor,
                                      const std::vector<int> &l_coor,
                                      const double u[],
                                      double M[])
{
    int max_ao_geo_order = 5; // FIXME hardcoded
    int buffer_len = balboa_get_buffer_len(balboa_context, max_ao_geo_order, AO_BLOCK_LENGTH);
    double *k_ao_compressed = new double[buffer_len];
    int *k_ao_compressed_index = new int[buffer_len];
    int k_ao_compressed_num;
    double *l_ao_compressed = new double[buffer_len];
    int *l_ao_compressed_index = new int[buffer_len];
    int l_ao_compressed_num;
    int num_aos = balboa_get_num_aos(balboa_context);
    int *ao_centers = new int[num_aos];
    for (int i = 0; i < num_aos; i++)
    {
        ao_centers[i] = balboa_get_ao_center(balboa_context, i);
    }
    int slice_offsets[4];
    compute_slice_offsets(k_coor, slice_offsets);
    compress(use_gradient,
             k_ao_compressed_num,
             k_ao_compressed_index,
             k_ao_compressed,
             num_aos,
             ao,
             ao_centers,
             k_coor,
             slice_offsets);
    compute_slice_offsets(l_coor, slice_offsets);
    compress(use_gradient,
             l_ao_compressed_num,
             l_ao_compressed_index,
             l_ao_compressed,
             num_aos,
             ao,
             ao_centers,
             l_coor,
             slice_offsets);

    double prefactors[5] = {f, f, f, f, 0.5 * f};

    // distribute over XC potential derivative matrix
    distribute_matrix(mat_dim,
                      use_gradient,
                      use_tau,
                      prefactors,
                      u,
                      M,
                      k_ao_compressed_num,
                      k_ao_compressed_index,
                      k_ao_compressed,
                      l_ao_compressed_num,
                      l_ao_compressed_index,
                      l_ao_compressed);
    if (use_gradient)
    {
        prefactors[0] = 0.0;
        distribute_matrix(mat_dim,
                          use_gradient,
                          false,
                          prefactors,
                          u,
                          M,
                          l_ao_compressed_num,
                          l_ao_compressed_index,
                          l_ao_compressed,
                          k_ao_compressed_num,
                          k_ao_compressed_index,
                          k_ao_compressed);
    }
    delete[] k_ao_compressed;
    delete[] k_ao_compressed_index;
    delete[] l_ao_compressed;
    delete[] l_ao_compressed_index;
    delete[] ao_centers;
}

int AOBatch::set_basis(const int basis_type,
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

int AOBatch::get_num_aos() { return balboa_get_num_aos(balboa_context); }


void AOBatch::compute_slice_offsets(const std::vector<int> &coor,
                                    int off[])
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
