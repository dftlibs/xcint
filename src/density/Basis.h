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
                  const double in_center_xyz[],
                  const int    in_center_element[],
                  const int    in_num_shells,
                  const int    in_shell_center[],
                  const int    in_l_quantum_num[],
                  const int    in_shell_num_primitives[],
                  const double in_primitive_exp[],
                  const double in_contraction_coef[]);
        int  get_num_centers() const;
        int  get_num_ao_slices() const;
        int  get_num_ao() const;
        int  get_num_ao_cartesian() const;
        int  get_ao_center(const int i) const;
        int  get_geo_off(const int i,
                         const int j,
                         const int k) const;
        void report();
        void set_geo_off(const int geo_diff_order);

        int     num_centers; // FIXME
        int     num_shells; // FIXME
        int    *l_quantum_num; // FIXME
        int    *shell_num_primitives; // FIXME
        double *primitive_exp; // FIXME
        double *center_xyz; // FIXME
        int    *center_element; // FIXME
        int    *shell_center; // FIXME

    private:

        Basis(const Basis &rhs);            // not implemented
        Basis &operator=(const Basis &rhs); // not implemented

        void    nullify();

        double *shell_center_xyz;
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
        double *contraction_coef;
        bool    is_initialized;
        bool    is_synced;
};

#endif // Basis_h_
