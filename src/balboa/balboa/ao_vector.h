#pragma once

// vec r = vec p * vec a
void get_pa(const int num_points,
            const double *__restrict__ a,
            const double *__restrict__ p,
            double *__restrict__ r);

// vec r = vec p * vec a
void get_pa_block(const double *__restrict__ a,
                  const double *__restrict__ p,
                  double *__restrict__ r);

// vec r = vec p * vec a + vec b
void get_pa_plus_b(const int num_points,
                   const double *__restrict__ a,
                   const double *__restrict__ p,
                   const double *__restrict__ b,
                   double *__restrict__ r);

// vec r = vec p * vec a + vec b
void get_pa_plus_b_block(const double *__restrict__ a,
                         const double *__restrict__ p,
                         const double *__restrict__ b,
                         double *__restrict__ r);

// vec r = vec p * vec a + s * vec b
void get_pa_plus_sb(const int num_points,
                    const double *__restrict__ a,
                    const double *__restrict__ p,
                    const double s,
                    const double *__restrict__ b,
                    double *__restrict__ r);

// vec r = vec p * vec a + s * vec b
void get_pa_plus_sb_block(const double *__restrict__ a,
                          const double *__restrict__ p,
                          const double s,
                          const double *__restrict__ b,
                          double *__restrict__ r);

void get_exp(const int num_points,
             const double *__restrict__ p2,
             const double c,
             const double a,
             double *__restrict__ s);

void get_exp_block(const double *__restrict__ p2,
                   const double c,
                   const double a,
                   double *__restrict__ s);

void vec_daxpy(const int num_points,
               const double s,
               const double *__restrict__ a,
               double *__restrict__ r);

void vec_daxpy_block(const double s,
                     const double *__restrict__ a,
                     double *__restrict__ r);

void get_p2(const int num_points,
            const double *__restrict__ shell_centers_coordinates,
            const double *__restrict__ x_coordinates_bohr,
            const double *__restrict__ y_coordinates_bohr,
            const double *__restrict__ z_coordinates_bohr,
            double *__restrict__ px,
            double *__restrict__ py,
            double *__restrict__ pz,
            double *__restrict__ p2);

void get_p2_block(const double *__restrict__ shell_centers_coordinates,
                  const double *__restrict__ x_coordinates_bohr,
                  const double *__restrict__ y_coordinates_bohr,
                  const double *__restrict__ z_coordinates_bohr,
                  double *__restrict__ px,
                  double *__restrict__ py,
                  double *__restrict__ pz,
                  double *__restrict__ p2);

bool calculate_chunk(const int num_points,
                     const double extent_squared,
                     const double *__restrict__ p2);

bool calculate_chunk_block(const double extent_squared,
                           const double *__restrict__ p2);
