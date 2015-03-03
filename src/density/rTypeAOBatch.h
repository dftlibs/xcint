#ifndef rTypeAOBatch_h_
#define rTypeAOBatch_h_

#include <stdlib.h>
#include <vector>

#include "rTypeBasis.h"

class rTypeAOBatch
{
    public:

        rTypeAOBatch();
        ~rTypeAOBatch();

        void get_ao(const rTypeBasis &basis,
                    const bool       use_gradient,
                    const int        max_ao_geo_order,
                    const double     p[]);

        void get_ao_shell(const int        ishell,
                          const rTypeBasis &basis,
                                double     ao_local[],
                          const int        max_ao_geo_order,
                          const double     p[]);

        void get_ao_shell(const int        ishell,
                          const rTypeBasis &basis,
                          const int        max_ao_geo_order,
                          const double     p[]);

        void distribute_matrix_undiff(const int    mat_dim,
                                      const bool   use_gradient,
                                      const bool   use_tau,
                                      const double prefactors[],
                                      const double u[],
                                            double fmat[]);

        void distribute_matrix(const int    mat_dim,
                               const bool   use_gradient,
                               const bool   use_tau,
                               const double prefactors[],
                               const double u[],
                                     double fmat[],
                               const int    k_aoc_num,
                               const int    k_aoc_index[],
                               const double k_aoc[],
                               const int    l_aoc_num,
                               const int    l_aoc_index[],
                               const double l_aoc[]);

        void get_density_undiff(const int    mat_dim,
                                const bool   use_gradient,
                                const bool   use_tau,
                                const double prefactors[],
                                      double density[],
                                const double dmat[],
                                const bool   dmat_is_symmetric,
                                const bool   kl_match);

        void get_density(const int    mat_dim,
                         const bool   use_gradient,
                         const bool   use_tau,
                         const double prefactors[],
                               double density[],
                         const double dmat[],
                         const bool   dmat_is_symmetric,
                         const bool   kl_match,
                         const int    k_aoc_num,
                         const int    k_aoc_index[],
                         const double k_aoc[],
                         const int    l_aoc_num,
                         const int    l_aoc_index[],
                         const double l_aoc[]);

        void get_dens_geo_derv(const rTypeBasis       &basis,
                               const int              mat_dim,
                               const bool             use_gradient,
                               const bool             use_tau,
                               const std::vector<int> &coor,
                                     double           density[],
                               const bool             get_dens,
                                     double           mat[]);

    private:

        rTypeAOBatch(const rTypeAOBatch &rhs);            // not implemented
        rTypeAOBatch &operator=(const rTypeAOBatch &rhs); // not implemented

        bool is_same_center(const int c,
                            const std::vector<int> &carray) const;

        void diff_wrt_center_tuple(const rTypeBasis       &basis,
                                   const int              mat_dim,
                                   const bool             use_gradient,
                                   const bool             use_tau,
                                   const double           f,
                                   const std::vector<int> &k_coor,
                                   const std::vector<int> &l_coor,
                                   const bool             get_dens,
                                         double           u[],
                                         double           M[]);

        void compress(const rTypeBasis       &basis,
                      const bool             use_gradient,
                            int              &aoc_num,
                            int              *(&aoc_index),
                            double           *(&aoc),
                      const std::vector<int> &coor);

        void transform_basis(const rTypeBasis &basis) const;

        void nullify();

        int     ao_length;
        double *ao;

        double *ao_compressed;
        int     ao_compressed_num;
        int    *ao_compressed_index;

        double *k_ao_compressed;
        int     k_ao_compressed_num;
        int    *k_ao_compressed_index;

        double *l_ao_compressed;
        int     l_ao_compressed_num;
        int    *l_ao_compressed_index;
};

#endif // rTypeAOBatch_h_
