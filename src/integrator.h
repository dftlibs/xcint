#pragma once

#include "Functional.h"
#include "balboa.h"
#include "xcint.h"

class XCint
{
  public:
    XCint();
    ~XCint();

    int set_basis(const int basis_type,
                  const int num_centers,
                  const double center_coordinates[],
                  const int num_shells,
                  const int shell_centers[],
                  const int shell_l_quantum_numbers[],
                  const int shell_num_primitives[],
                  const double primitive_exponents[],
                  const double contraction_coefficients[]);

    int set_functional(const char *line);

    int integrate(const xcint_mode_t mode,
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
                  double *num_electrons);
    //   double *num_electrons) const;

  private:
    XCint(const XCint &rhs);            // not implemented
    XCint &operator=(const XCint &rhs); // not implemented

    Functional fun;
    balboa_context_t *balboa_context;
    int ao_length;
    double *ao;

    void nullify();

    void distribute_matrix2(const int block_length,
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
                            const double grid_w[]);
    //            const double grid_w[]) const;

    void integrate_batch(const double dmat[],
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
                         //     const double grid_w[]) const;
                         const double grid_w[]);

    void compute_slice_offsets(const std::vector<int> &coor, int off[]);
};
