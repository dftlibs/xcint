#ifndef XCint_h_
#define XCint_h_

#include "Basis.h"
#include "AOBatch.h"
#include "Functional.h"

typedef int (*print_function)(const char* line);

class XCint
{
    public:

        XCint();
        ~XCint();

        int set_basis(const int    basis_type,
                      const int    num_centers,
                      const double center_coordinates[],
                      const int    center_elements[],
                      const int    num_shells,
                      const int    shell_centers[],
                      const int    shell_l_quantum_numbers[],
                      const int    shell_num_primitives[],
                      const double primitive_exponents[],
                      const double contraction_coefficients[]);

        int set_functional(const char *line);

        int integrate(const int    mode,
                      const int    num_points,
                      const double grid_pw[],
                      const int    num_pert,
                      const int    pert[],
                      const int    comp[],
                      const int    num_dmat,
                      const int    dmat_to_pert[],
                      const int    dmat_to_comp[],
                      const double dmat[],
                      const int    get_xc_energy,
                            double *xc_energy,
                      const int    get_xc_mat,
                            double xc_mat[],
                            double *num_electrons) const;

    private:

        XCint(const XCint &rhs);            // not implemented
        XCint &operator=(const XCint &rhs); // not implemented

        Functional fun;
        Basis basis;

        void nullify();

        void distribute_matrix(const int              block_length,
                               const int              num_variables,
                               const int              num_pert,
                               const int              mat_dim,
                               const double           prefactors[],
                               const int              w_off,
                               const bool             n_is_used[],
                               const double           n[],
                                     double           u[],
                                     double           xc_mat[],
                                     double           &xc_energy,
                               const std::vector<int> coor,
                                     AOBatch     &batch,
                               const double           grid_pw[]) const;

        void integrate_batch(const double dmat[],
                             const int    get_xc_energy,
                                   double &xc_energy,
                             const int    get_xc_mat,
                                   double xc_mat[],
                                   double &num_electrons,
                             const int    geo_coor[],
                             const bool   use_dmat[],
                             const int    ipoint,
                             const int    geo_derv_order,
                             const int    max_ao_order_g,
                             const int    block_length,
                             const int    num_variables,
                             const int    num_pert,
                             const int    num_fields,
                             const int    mat_dim,
                             const bool   get_gradient,
                             const bool   get_tau,
                             const int    dmat_index[],
                             const double grid_pw[]) const;

//      double time_total;
//      double time_ao;
//      double time_fun_derv;
//      double time_densities;
//      double time_matrix_distribution;

        void reset_time();
};

#endif // XCint_h_
