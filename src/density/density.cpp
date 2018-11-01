#include "blas_interface.h"
#include "compress.h"
#include "density.h"

// size_t
#include <cstddef>

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cassert>
#include <iostream>

void distribute_matrix(const int mat_dim,
                       const int block_length,
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

    double *W = new double[k_aoc_num * block_length];

    std::fill(&W[0], &W[block_length * k_aoc_num], 0.0);

    int kc, lc, iboff;

    int num_slices;
    (use_gradient) ? (num_slices = 4) : (num_slices = 1);

    for (int islice = 0; islice < num_slices; islice++)
    {
        if (std::abs(prefactors[islice]) > 0.0)
        {
            for (int k = 0; k < k_aoc_num; k++)
            {
                iboff = k * block_length;
                for (int ib = 0; ib < block_length; ib++)
                {
                    W[iboff + ib] +=
                        prefactors[islice] * u[islice * block_length + ib] *
                        k_aoc[islice * block_length * mat_dim + iboff + ib];
                }
            }
        }
    }

    double *F = new double[k_aoc_num * l_aoc_num];

    // we compute F(k, l) += W(k, b) AO_l(l, b)^T
    // we transpose W instead of AO_l because we call fortran blas
    char ta = 't';
    char tb = 'n';
    int im = l_aoc_num;
    int in = k_aoc_num;
    int ik = block_length;
    int lda = ik;
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
               W,
               &ldb,
               &beta,
               F,
               &ldc);

    if (use_tau)
    {
        if (std::abs(prefactors[4]) > 0.0)
        {
            for (int ixyz = 0; ixyz < 3; ixyz++)
            {
                for (int k = 0; k < k_aoc_num; k++)
                {
                    iboff = k * block_length;
                    for (int ib = 0; ib < block_length; ib++)
                    {
                        W[iboff + ib] =
                            u[4 * block_length + ib] *
                            k_aoc[(ixyz + 1) * block_length * mat_dim + iboff +
                                  ib];
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
                           &l_aoc[(ixyz + 1) * block_length * mat_dim],
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

void get_density(const int mat_dim,
                 const int block_length,
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
    double *X = new double[k_aoc_num * block_length];

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
        int im = block_length;
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
        int im = block_length;
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

//  it is not ok to zero out here since geometric
//  derivatives are accumulated from several contributions
//  std::fill(&density[0], &density[num_slices * block_length], 0.0);

    // assemble density and possibly gradient
    for (int islice = 0; islice < num_slices; islice++)
    {
        if (std::abs(prefactors[islice]) > 0.0)
        {
            for (int k = 0; k < k_aoc_num; k++)
            {
                iboff = k * block_length;
                for (int ib = 0; ib < block_length; ib++)
                {
                    density[islice * block_length + ib] +=
                        prefactors[islice] * X[iboff + ib] *
                        k_aoc[islice * block_length * mat_dim + iboff + ib];
                }
            }
        }
    }

    if (use_tau)
    {
        if (std::abs(prefactors[4]) > 0.0)
        {
            for (int ixyz = 0; ixyz < 3; ixyz++)
            {
                if (dmat_is_symmetric)
                {
                    char si = 'r';
                    char up = 'u';
                    int im = block_length;
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
                               &l_aoc[(ixyz + 1) * block_length * mat_dim],
                               &ldb,
                               &beta,
                               X,
                               &ldc);
                }
                else
                {
                    char ta = 'n';
                    char tb = 'n';
                    int im = block_length;
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
                               &l_aoc[(ixyz + 1) * block_length * mat_dim],
                               &lda,
                               D,
                               &ldb,
                               &beta,
                               X,
                               &ldc);
                }

                for (int k = 0; k < k_aoc_num; k++)
                {
                    iboff = k * block_length;
                    for (int ib = 0; ib < block_length; ib++)
                    {
                        density[4 * block_length + ib] +=
                            prefactors[4] * X[iboff + ib] *
                            k_aoc[(ixyz + 1) * block_length * mat_dim + iboff +
                                  ib];
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

void get_dens_geo_derv(const int mat_dim,
                       const int num_aos,
                       const int block_length,
                       const int buffer_len,
                       const double ao[],
                       const int ao_centers[],
                       const bool use_gradient,
                       const bool use_tau,
                       const std::vector<int> &coor,
                       std::function<int(int, int, int)> get_geo_offset,
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
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 2:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 3:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[1]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 4:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        k_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[1]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[3]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_u_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        break;
    default:
        std::cout << "ERROR: get_dens_geo_derv coor.size() too high\n";
        exit(1);
        break;
    }
}

void get_mat_geo_derv(const int mat_dim,
                      const int num_aos,
                      const int block_length,
                      const int buffer_len,
                      const double ao[],
                      const int ao_centers[],
                      const bool use_gradient,
                      const bool use_tau,
                      const std::vector<int> &coor,
                      std::function<int(int, int, int)> get_geo_offset,
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
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 2:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 3:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[1]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        break;
    case 4:
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        k_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        l_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[1]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[1]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[2]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[3]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        k_coor.push_back(coor[0]);
        k_coor.push_back(coor[3]);
        l_coor.push_back(coor[1]);
        l_coor.push_back(coor[2]);
        diff_M_wrt_center_tuple(mat_dim,
                                num_aos,
                                block_length,
                                buffer_len,
                                ao,
                                ao_centers,
                                use_gradient,
                                use_tau,
                                f,
                                get_geo_offset,
                                k_coor,
                                l_coor,
                                density,
                                mat);
        k_coor.clear();
        l_coor.clear();
        break;
    default:
        std::cout << "ERROR: get_dens_geo_derv coor.size() too high\n";
        exit(1);
        break;
    }
}

void diff_u_wrt_center_tuple(const int mat_dim,
                             const int num_aos,
                             const int block_length,
                             const int buffer_len,
                             const double ao[],
                             const int ao_centers[],
                             const bool use_gradient,
                             const bool use_tau,
                             const double f,
                             std::function<int(int, int, int)> get_geo_offset,
                             const std::vector<int> &k_coor,
                             const std::vector<int> &l_coor,
                             double u[],
                             const double M[])
{
    double *k_ao_compressed = new double[buffer_len];
    int *k_ao_compressed_index = new int[buffer_len];
    int k_ao_compressed_num;
    double *l_ao_compressed = new double[buffer_len];
    int *l_ao_compressed_index = new int[buffer_len];
    int l_ao_compressed_num;
    int slice_offsets[4];
    compute_slice_offsets(get_geo_offset, k_coor, slice_offsets);
    compress(use_gradient,
             block_length,
             k_ao_compressed_num,
             k_ao_compressed_index,
             k_ao_compressed,
             num_aos,
             ao,
             ao_centers,
             k_coor,
             slice_offsets);
    compute_slice_offsets(get_geo_offset, l_coor, slice_offsets);
    compress(use_gradient,
             block_length,
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
                block_length,
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
                    block_length,
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
}

void diff_M_wrt_center_tuple(const int mat_dim,
                             const int num_aos,
                             const int block_length,
                             const int buffer_len,
                             const double ao[],
                             const int ao_centers[],
                             const bool use_gradient,
                             const bool use_tau,
                             const double f,
                             std::function<int(int, int, int)> get_geo_offset,
                             const std::vector<int> &k_coor,
                             const std::vector<int> &l_coor,
                             const double u[],
                             double M[])
{
    double *k_ao_compressed = new double[buffer_len];
    int *k_ao_compressed_index = new int[buffer_len];
    int k_ao_compressed_num;
    double *l_ao_compressed = new double[buffer_len];
    int *l_ao_compressed_index = new int[buffer_len];
    int l_ao_compressed_num;
    int slice_offsets[4];
    compute_slice_offsets(get_geo_offset, k_coor, slice_offsets);
    compress(use_gradient,
             block_length,
             k_ao_compressed_num,
             k_ao_compressed_index,
             k_ao_compressed,
             num_aos,
             ao,
             ao_centers,
             k_coor,
             slice_offsets);
    compute_slice_offsets(get_geo_offset, l_coor, slice_offsets);
    compress(use_gradient,
             block_length,
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
                      block_length,
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
                          block_length,
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
}

void compute_slice_offsets(std::function<int(int, int, int)> get_geo_offset,
                           const std::vector<int> &coor,
                           int off[])
{
    int kp[3] = {0, 0, 0};
    for (size_t j = 0; j < coor.size(); j++)
    {
        kp[(coor[j] - 1) % 3]++;
    }
    off[0] = get_geo_offset(kp[0], kp[1], kp[2]);
    off[1] = get_geo_offset(kp[0] + 1, kp[1], kp[2]);
    off[2] = get_geo_offset(kp[0], kp[1] + 1, kp[2]);
    off[3] = get_geo_offset(kp[0], kp[1], kp[2] + 1);
}
