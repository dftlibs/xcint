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
                                   const char      *line);

typedef enum
{
    XCINT_BASIS_SPHERICAL,
    XCINT_BASIS_CARTESIAN
} xcint_basis_t;

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
                              const double        contraction_coefficients[]);

typedef enum
{
    XCINT_MODE_RKS,
    XCINT_MODE_UKS
} xcint_mode_t;

typedef enum
{
    XCINT_PERT_EL,
    XCINT_PERT_GEO,
    XCINT_PERT_MAG_CGO,
    XCINT_PERT_MAG_LAO
} xcint_perturbation_t;

XCINT_API int xcint_integrate(const xcint_context_t      *context,
                              const xcint_mode_t         mode,
                              const int                  num_points,
                              const double               grid[],
                              const int                  num_perturbations,
                              const xcint_perturbation_t perturbations[],
                              const int                  components[],
                              const int                  num_dmat,
                              const int                  dmat_to_perturbations[],
                              const double               dmat[],
                              const int                  get_exc,
                                    double               *exc,
                              const int                  get_vxc,
                                    double               vxc[],
                                    double               *num_electrons);

#ifdef __cplusplus
}
#endif

#endif /* XCINT_H_INCLUDED */
