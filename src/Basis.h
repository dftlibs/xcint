#ifndef Basis_h_
#define Basis_h_

class Basis
{
    friend class AOBatch;

    public:

        Basis();
        ~Basis();

        void init(const int    in_basis_type,
                  const int    in_num_centers,
                  const double in_center_coordinates[],
                  const int    in_num_shells,
                  const int    in_shell_centers[],
                  const int    in_shell_l_quantum_numbers[],
                  const int    in_shell_num_primitives[],
                  const double in_primitive_exponents[],
                  const double in_contraction_coefficients[]);
        int  get_num_centers() const;
        int  get_num_ao_slices() const;
        int  get_num_ao() const;
        int  get_num_ao_cartesian() const;
        int  get_ao_center(const int i) const;
        int  get_geo_off(const int i,
                         const int j,
                         const int k) const;
        void set_geo_off(const int geo_diff_order);

        int     num_centers; // FIXME
        int     num_shells; // FIXME
        int    *shell_l_quantum_numbers; // FIXME
        int    *shell_num_primitives; // FIXME
        double *primitive_exponents; // FIXME
        double *center_coordinates; // FIXME
        int    *shell_centers; // FIXME

    private:

        Basis(const Basis &rhs);            // not implemented
        Basis &operator=(const Basis &rhs); // not implemented

        void    nullify();
        void    deallocate();

        double *shell_centers_coordinates;
        double *shell_extent_squared;
        int    *cartesian_deg;
        int    *shell_off;
        int    *spherical_deg;
        bool    is_spherical;
        int     num_ao;
        int     num_ao_cartesian;
        int     num_ao_spherical;
        int     num_ao_slices;
        int    *ao_center;
        int     geo_diff_order;
        int    *geo_off;
        double *contraction_coefficients;
        int     is_initialized;
};

#endif // Basis_h_
