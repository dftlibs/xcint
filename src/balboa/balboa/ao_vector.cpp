#include <math.h>

#include "ao_vector.h"
#include "balboa_parameters.h"

// vec r = vec p * vec a
void get_pa(const int num_points,
            const double *__restrict__ a,
            const double *__restrict__ p,
            double *__restrict__ r)
{
    for (int k = 0; k < num_points; k++)
    {
        r[k] = p[k] * a[k];
    }
}

// vec r = vec p * vec a
void get_pa_block(const double *__restrict__ a,
                  const double *__restrict__ p,
                  double *__restrict__ r)
{
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        r[k] = p[k] * a[k];
    }
}

// vec r = vec p * vec a + vec b
void get_pa_plus_b(const int num_points,
                   const double *__restrict__ a,
                   const double *__restrict__ p,
                   const double *__restrict__ b,
                   double *__restrict__ r)
{
    for (int k = 0; k < num_points; k++)
    {
        r[k] = p[k] * a[k] + b[k];
    }
}

// vec r = vec p * vec a + vec b
void get_pa_plus_b_block(const double *__restrict__ a,
                         const double *__restrict__ p,
                         const double *__restrict__ b,
                         double *__restrict__ r)
{
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        r[k] = p[k] * a[k] + b[k];
    }
}

// vec r = vec p * vec a + s * vec b
void get_pa_plus_sb(const int num_points,
                    const double *__restrict__ a,
                    const double *__restrict__ p,
                    const double s,
                    const double *__restrict__ b,
                    double *__restrict__ r)
{
    for (int k = 0; k < num_points; k++)
    {
        r[k] = p[k] * a[k] + s * b[k];
    }
}

// vec r = vec p * vec a + s * vec b
void get_pa_plus_sb_block(const double *__restrict__ a,
                          const double *__restrict__ p,
                          const double s,
                          const double *__restrict__ b,
                          double *__restrict__ r)
{
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        r[k] = p[k] * a[k] + s * b[k];
    }
}

void get_exp(const int num_points,
             const double *__restrict__ p2,
             const double c,
             const double a,
             double *__restrict__ s)
{
    double b1, b2;

    for (int k = 0; k < num_points; k++)
    {
        b1 = a * p2[k];
        b2 = exp(b1);
        s[k] = c * b2;
    }
}

void get_exp_block(const double *__restrict__ p2,
                   const double c,
                   const double a,
                   double *__restrict__ s)
{
    double b1, b2;

    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        b1 = a * p2[k];
        b2 = exp(b1);
        s[k] = c * b2;
    }
}

void vec_daxpy(const int num_points,
               const double s,
               const double *__restrict__ a,
               double *__restrict__ r)
{
    for (int k = 0; k < num_points; k++)
    {
        r[k] += s * a[k];
    }
}

void vec_daxpy_block(const double s,
                     const double *__restrict__ a,
                     double *__restrict__ r)
{
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        r[k] += s * a[k];
    }
}

void get_p2(const int num_points,
            const double *__restrict__ shell_centers_coordinates,
            const double *__restrict__ x_coordinates_bohr,
            const double *__restrict__ y_coordinates_bohr,
            const double *__restrict__ z_coordinates_bohr,
            double *__restrict__ px,
            double *__restrict__ py,
            double *__restrict__ pz,
            double *__restrict__ p2)
{
    for (int k = 0; k < num_points; k++)
    {
        px[k] = x_coordinates_bohr[k] - shell_centers_coordinates[0];
        py[k] = y_coordinates_bohr[k] - shell_centers_coordinates[1];
        pz[k] = z_coordinates_bohr[k] - shell_centers_coordinates[2];
        p2[k] = px[k] * px[k] + py[k] * py[k] + pz[k] * pz[k];
    }
}

void get_p2_block(const double *__restrict__ shell_centers_coordinates,
                  const double *__restrict__ x_coordinates_bohr,
                  const double *__restrict__ y_coordinates_bohr,
                  const double *__restrict__ z_coordinates_bohr,
                  double *__restrict__ px,
                  double *__restrict__ py,
                  double *__restrict__ pz,
                  double *__restrict__ p2)
{
    for (int k = 0; k < AO_CHUNK_LENGTH; k++)
    {
        px[k] = x_coordinates_bohr[k] - shell_centers_coordinates[0];
        py[k] = y_coordinates_bohr[k] - shell_centers_coordinates[1];
        pz[k] = z_coordinates_bohr[k] - shell_centers_coordinates[2];
        p2[k] = px[k] * px[k] + py[k] * py[k] + pz[k] * pz[k];
    }
}

bool calculate_chunk(const int num_points,
                     const double extent_squared,
                     const double *__restrict__ p2)
{
    for (int k = 0; k < num_points; k++)
    {
        if (p2[k] < extent_squared)
        {
            return true;
        }
    }
    return false;
}

bool calculate_chunk_block(const double extent_squared,
                           const double *__restrict__ p2)
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
