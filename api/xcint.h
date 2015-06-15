#ifndef XCINT_H_INCLUDED
#define XCINT_H_INCLUDED

#ifndef XCINT_API
#  ifdef _WIN32
#     if defined(XCINT_BUILD_SHARED) /* build dll */
#         define XCINT_API __declspec(dllexport)
#     elif !defined(XCINT_BUILD_STATIC) /* use dll */
#         define XCINT_API __declspec(dllimport)
#     else /* static library */
#         define XCINT_API
#     endif
#  else
#     if __GNUC__ >= 4
#         define XCINT_API __attribute__((visibility("default")))
#     else
#         define XCINT_API
#     endif
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct xcint_context_s;
typedef struct xcint_context_s xcint_context_t;

XCINT_API xcint_context_t *xcint_new();
XCINT_API void xcint_free(xcint_context_t *context);

XCINT_API int xcint_set_functional(xcint_context_t *context,
                                   const char   *line,
                                         double &hfx,
                                         double &mu,
                                         double &beta);

XCINT_API void xcint_set_basis(xcint_context_t *context,
                               const int    basis_type,
                               const int    num_centers,
                               const double center_xyz[],
                               const int    center_element[],
                               const int    num_shells,
                               const int    shell_center[],
                               const int    l_quantum_num[],
                               const int    shell_num_primitives[],
                               const double primitive_exp[],
                               const double contraction_coef[]);

// FIXME should be const
XCINT_API void xcint_integrate(xcint_context_t *context,
                               const int    mode,
                               const int    num_points,
                               const double grid_pw[],
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

#ifdef __cplusplus
}
#endif

#endif
