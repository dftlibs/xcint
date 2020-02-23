#pragma once

class Main
{
  public:
    Main();
    ~Main();

    int set_basis(const int basis_type,
                  const int num_centers,
                  const double center_coordinates_bohr[],
                  const int num_shells,
                  const int shell_centers[],
                  const int shell_l_quantum_numbers[],
                  const int shell_num_primitives[],
                  const double primitive_exponents[],
                  const double contraction_coefficients[]);

    int get_buffer_len(const int max_geo_order, const int num_points) const;
    int get_ao_center(const int i) const;
    int get_geo_offset(const int i, const int j, const int k) const;

    int get_num_aos() const;

    // buffer is not zeroed out inside get_ao
    int get_ao(const int max_geo_order,
               const int num_points,
               const double x_coordinates_bohr[],
               const double y_coordinates_bohr[],
               const double z_coordinates_bohr[],
               double buffer[]) const;

  private:
    Main(const Main &rhs);            // not implemented
    Main &operator=(const Main &rhs); // not implemented

    void nullify();
    void deallocate();

    void transform_basis() const;

    int num_centers;
    int num_shells;
    int *shell_l_quantum_numbers;
    int *shell_num_primitives;
    double *primitive_exponents;
    double *center_coordinates_bohr;
    int *shell_centers;

    double *shell_centers_coordinates;
    double *shell_extent_squared;
    int *cartesian_deg;
    int *shell_off;
    int *spherical_deg;
    bool is_spherical;
    int num_ao;
    int num_ao_cartesian;
    int num_ao_spherical;
    int num_ao_slices;
    int *ao_center;
    double *contraction_coefficients;
    int is_initialized;
    int *geo_offset;
    int geo_offset_size;
};
