#ifndef AOBatch_h_
#define AOBatch_h_

#include <stdlib.h>
#include <vector>

#include "Basis.h"
#include "balboa.h"

class AOBatch
{
  public:
    AOBatch();
    ~AOBatch();

    int set_basis(const int basis_type,
                  const int num_centers,
                  const double center_coordinates_bohr[],
                  const int num_shells,
                  const int shell_centers[],
                  const int shell_l_quantum_numbers[],
                  const int shell_num_primitives[],
                  const double primitive_exponents[],
                  const double contraction_coefficients[]);

    void get_ao(const Basis &basis,
                const bool use_gradient,
                const int max_ao_geo_order,
                const int block_length,
                const double p[]);

    int get_num_aos();

    void distribute_matrix_undiff(const int mat_dim,
                                  const bool use_gradient,
                                  const bool use_tau,
                                  const double prefactors[],
                                  const double u[],
                                  double fmat[]);

    void distribute_matrix(const int mat_dim,
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
                           const double l_aoc[]);

    void get_density_undiff(const int mat_dim,
                            const bool use_gradient,
                            const bool use_tau,
                            const double prefactors[],
                            double density[],
                            const double dmat[],
                            const bool dmat_is_symmetric,
                            const bool kl_match);

    void get_density(const int mat_dim,
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
                     const double l_aoc[]);

    void get_mat_geo_derv(const Basis &basis,
                          const int mat_dim,
                          const bool use_gradient,
                          const bool use_tau,
                          const std::vector<int> &coor,
                          const double density[],
                          double mat[]);

    void get_dens_geo_derv(const Basis &basis,
                           const int mat_dim,
                           const bool use_gradient,
                           const bool use_tau,
                           const std::vector<int> &coor,
                           double density[],
                           const double mat[]);

  private:
    AOBatch(const AOBatch &rhs);            // not implemented
    AOBatch &operator=(const AOBatch &rhs); // not implemented

    balboa_context_t *balboa_context;

    bool is_same_center(const int c, const std::vector<int> &carray);

    void diff_M_wrt_center_tuple(const Basis &basis,
                                 const int mat_dim,
                                 const bool use_gradient,
                                 const bool use_tau,
                                 const double f,
                                 const std::vector<int> &k_coor,
                                 const std::vector<int> &l_coor,
                                 const double u[],
                                 double M[]);

    void diff_u_wrt_center_tuple(const Basis &basis,
                                 const int mat_dim,
                                 const bool use_gradient,
                                 const bool use_tau,
                                 const double f,
                                 const std::vector<int> &k_coor,
                                 const std::vector<int> &l_coor,
                                 double u[],
                                 const double M[]);

    void compress(const Basis &basis,
                  const bool use_gradient,
                  int &aoc_num,
                  int *(&aoc_index),
                  double *(&aoc),
                  const std::vector<int> &coor);

    void transform_basis(const Basis &basis);

    void nullify();

    int ao_length;
    double *ao;

    double *ao_compressed;
    int ao_compressed_num;
    int *ao_compressed_index;

    double *k_ao_compressed;
    int k_ao_compressed_num;
    int *k_ao_compressed_index;

    double *l_ao_compressed;
    int l_ao_compressed_num;
    int *l_ao_compressed_index;
};

#endif // AOBatch_h_
