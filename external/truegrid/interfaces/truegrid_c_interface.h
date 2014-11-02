#ifndef truegrid_c_interface_h_
#define truegrid_c_interface_h_

extern "C"
{
void truegrid_set_grid_parameters(const double radial_precision,
                                  const int    angular_min,
                                  const int    angular_max);

void truegrid_generate(const int    verbosity,
                       const int    num_centers,
                       const double center_xyz[],
                       const int    center_element[],
                       const int    num_shells,
                       const int    shell_center[],
                       const int    l_quantum_num[],
                       const int    shell_num_primitives[],
                       const double primitive_exp[]);

void truegrid_read();

void truegrid_stretch();

int truegrid_get_num_points();

double *truegrid_get_grid_w();

double *truegrid_get_grid_p();
}

#endif // truegrid_c_interface_h_
