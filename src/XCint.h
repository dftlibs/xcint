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

        void set_basis(const int    basis_type,
                       const int    num_centers,
                       const double center_xyz[],
                       const int    center_element[],
                       const int    num_shells,
                       const int    shell_center[],
                       const int    l_quantum_num[],
                       const int    shell_num_primitives[],
                       const double primitive_exp[],
                       const double contraction_coef[]);

        void print_splash();

        void set_verbosity(const int v);

        void set_stdout_function(print_function fun);
        void set_stderr_function(print_function fun);

        void set_functional(const char *line, double &hfx, double &mu, double &beta);

        void generate_grid(const double radial_precision,
                           const int    angular_min,
                           const int    angular_max,
                           const int    num_centers,
                           const double center_xyz[],
                           const int    center_element[],
                           const int    num_shells,
                           const int    shell_center[],
                           const int    l_quantum_num[],
                           const int    shell_num_primitives[],
                           const double primitive_exp[]);

        void integrate(const int    mode,
                       const int    num_pert,
                       const int    pert[],
                       const int    comp[],
                       const int    num_dmat,
                       const int    dmat_to_pert[],
                       const int    dmat_to_comp[],
                             double dmat[],
                       const int    get_xc_energy,
                             double &xc_energy,
                       const int    get_xc_mat,
                             double xc_mat[],
                             double &num_electrons);

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
                               const double           grid_w[]);

        void integrate_batch(      double dmat[],
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
                             const double grid_p[],
                             const double grid_w[]);

        double time_total;
        double time_ao;
        double time_fun_derv;
        double time_densities;
        double time_matrix_distribution;

        void reset_time();
        void print_timing();

        int verbosity;
};

#endif // XCint_h_
