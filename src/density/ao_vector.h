#ifndef ao_vector_h_
#define ao_vector_h_

// vec r = vec p * vec a
void get_pa(const double* __restrict__ a,
            const double* __restrict__ p,
                  double* __restrict__ r);

// vec r = vec p * vec a + vec b
void get_pa_plus_b(const double* __restrict__ a,
                   const double* __restrict__ p,
                   const double* __restrict__ b,
                         double* __restrict__ r);

// vec r = vec p * vec a + s * vec b
void get_pa_plus_sb(const double* __restrict__ a,
                    const double* __restrict__ p,
                    const double  s,
                    const double* __restrict__ b,
                          double* __restrict__ r);

void get_exp(const double* __restrict__ p2,
             const double  c,
             const double  a,
                   double* __restrict__ s);

void vec_daxpy(const double  s,
               const double* __restrict__ a,
                     double* __restrict__ r);

void get_p2(const double* __restrict__ shell_centers_coordinates,
            const double* __restrict__ p,
                  double* __restrict__ px,
                  double* __restrict__ py,
                  double* __restrict__ pz,
                  double* __restrict__ p2);

bool calculate_chunk(const double  extent_squared,
                     const double* __restrict__ p2);

#endif // ao_vector_h_
