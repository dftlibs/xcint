#ifndef xcint_c_api_h
#define xcint_c_api_h

#include "xcint_c_parameters.h"

extern "C"
{
    void xcint_print_splash();

    int xcint_set_functional(const char   *line,
                                   double &hfx,
                                   double &mu,
                                   double &beta);

    int xcint_generate_grid(const double radial_precision,
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

    void xcint_set_basis(const int    basis_type,
                         const int    num_centers,
                         const double center_xyz[],
                         const int    center_element[],
                         const int    num_shells,
                         const int    shell_center[],
                         const int    l_quantum_num[],
                         const int    shell_num_primitives[],
                         const double primitive_exp[],
                         const double contraction_coef[]);

    void xcint_integrate(const int    mode,
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

    void xcint_set_verbosity(const int v);

    void xcint_set_stdout_function(int (*fun)(const char* line));

    void xcint_set_stderr_function(int (*fun)(const char* line));
}

#endif // xcint_c_api_h
