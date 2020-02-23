#ifndef BALBOA_H_INCLUDED
#define BALBOA_H_INCLUDED

#ifndef BALBOA_API
#include "balboa_export.h"
#define BALBOA_API balboa_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct balboa_context_s;
typedef struct balboa_context_s balboa_context_t;

BALBOA_API
balboa_context_t *balboa_new_context();

BALBOA_API
void balboa_free_context(balboa_context_t *balboa_context);

BALBOA_API
int balboa_set_basis(balboa_context_t *balboa_context,
                     const int basis_type,
                     const int num_centers,
                     const double center_coordinates_bohr[],
                     const int num_shells,
                     const int shell_centers[],
                     const int shell_l_quantum_numbers[],
                     const int shell_num_primitives[],
                     const double primitive_exponents[],
                     const double contraction_coefficients[]);

BALBOA_API
int balboa_get_num_aos(const balboa_context_t *balboa_context);

BALBOA_API
int balboa_get_ao_center(const balboa_context_t *balboa_context, const int i);

BALBOA_API
int balboa_get_geo_offset(const balboa_context_t *balboa_context,
                          const int i,
                          const int j,
                          const int k);

BALBOA_API
int balboa_get_buffer_len(const balboa_context_t *balboa_context,
                          const int max_geo_order,
                          const int num_points);

BALBOA_API
int balboa_get_ao(const balboa_context_t *balboa_context,
                  const int max_geo_order,
                  const int num_points,
                  const double x_coordinates_bohr[],
                  const double y_coordinates_bohr[],
                  const double z_coordinates_bohr[],
                  double buffer[]);

#ifdef __cplusplus
}
#endif

#endif /* BALBOA_H_INCLUDED */
