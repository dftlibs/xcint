

#include <math.h>


#include "ao_vector.h"
#include "parameters.h"


// vec r = vec p * vec a
void get_pa(const double* __restrict__ a,
            const double* __restrict__ p,
                  double* __restrict__ r)
{
    #pragma ivdep
    #pragma vector aligned
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        r[k] = p[k]*a[k];
    }
}


// vec r = vec p * vec a + vec b
void get_pa_plus_b(const double* __restrict__ a,
                   const double* __restrict__ p,
                   const double* __restrict__ b,
                         double* __restrict__ r)
{
    #pragma ivdep
    #pragma vector aligned
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        r[k] = p[k]*a[k] + b[k];
    }
}


// vec r = vec p * vec a + s * vec b
void get_pa_plus_sb(const double* __restrict__ a,
                    const double* __restrict__ p,
                    const double  s,
                    const double* __restrict__ b,
                          double* __restrict__ r)
{
    #pragma ivdep
    #pragma vector aligned
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        r[k] = p[k]*a[k] + s*b[k];
    }
}


void get_exp(const double* __restrict__ p2,
             const double  c,
             const double  a,
                   double* __restrict__ s)
{
    double b1, b2;

    #pragma ivdep
    #pragma vector aligned
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        b1 = a*p2[k];
        b2 = exp(b1);
        s[k] = c*b2;
    }
}


void vec_daxpy(const double  s,
               const double* __restrict__ a,
                     double* __restrict__ r)
{
    #pragma ivdep
    #pragma vector aligned
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        r[k] += s*a[k];
    }
}


void get_p2(const double* __restrict__ shell_centers_coordinates,
            const double* __restrict__ pw,
                  double* __restrict__ px,
                  double* __restrict__ py,
                  double* __restrict__ pz,
                  double* __restrict__ p2)
{
    #pragma ivdep
    #pragma vector aligned
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        px[k] = pw[k*4    ] - shell_centers_coordinates[0];
        py[k] = pw[k*4 + 1] - shell_centers_coordinates[1];
        pz[k] = pw[k*4 + 2] - shell_centers_coordinates[2];
        p2[k] = px[k]*px[k] + py[k]*py[k] + pz[k]*pz[k];
    }
}


bool calculate_chunk(const double  extent_squared,
                     const double* __restrict__ p2)
{
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        if (p2[k] < extent_squared)
        {
            return true;
        }
    }
    return false;
}
