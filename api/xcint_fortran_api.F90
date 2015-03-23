module xcint_fortran_api

   use, intrinsic :: iso_c_binding
   implicit none

! definition of parameters
#include "xcint_fortran_parameters.h"

contains

   integer(c_int) function xcint_set_functional(line, &
                                                hfx,  &
                                                mu,   &
                                                beta) &
                           bind (c)
   ! input
      ! functional line (c_null_char-terminated)
      character(c_char), intent(in) :: line
   ! output
      ! orbital exchange factor ("HF exchange", "exact exchange")
      real(c_double) :: hfx
      ! cam mu
      real(c_double) :: mu
      ! cam beta
      real(c_double) :: beta
   ! returns
   !    0 upon success
   end function

   integer(c_int) function xcint_generate_grid(radial_precision,     &
                                               angular_min,          &
                                               angular_max,          &
                                               num_centers,          &
                                               center_xyz,           &
                                               center_element,       &
                                               num_shells,           &
                                               shell_center,         &
                                               l_quantum_num,        &
                                               shell_num_primitives, &
                                               primitive_exp)        &
                           bind (c)
   ! input
      ! desired radial precision
      real(c_double), intent(in), value :: radial_precision
      ! min number of angular points (pruning)
      integer(c_int), intent(in), value :: angular_min
      ! max number of angular points
      integer(c_int), intent(in), value :: angular_max
      ! number of centers (atoms)
      integer(c_int), value :: num_centers
      ! center coordinates, dimension is 3*num_centers
      real(c_double)        :: center_xyz(*)
      ! center element (1 for H, 2 for He, ...), dimension is num_centers
      integer(c_int)        :: center_element(*)
      ! number of shells
      integer(c_int), value :: num_shells
      ! shell center, dimension is num_shells
      integer(c_int)        :: shell_center(*)
      ! L quantum number for each shell, dimension is num_shells
      integer(c_int)        :: l_quantum_num(*)
      ! number of primitives for each shell, dimension is num_shells
      integer(c_int)        :: shell_num_primitives(*)
      ! primitive exponents, dimension is sum(shell_num_primitives)
      real(c_double)        :: primitive_exp(*)
   ! returns
   !    0 upon success
   end function

   subroutine xcint_set_basis(basis_type,           &
                              num_centers,          &
                              center_xyz,           &
                              center_element,       &
                              num_shells,           &
                              shell_center,         &
                              l_quantum_num,        &
                              shell_num_primitives, &
                              primitive_exp,        &
                              contraction_coef)     &
              bind (c)
   ! input
      ! basis set type (XCINT_BASIS_SPHERICAL or XCINT_BASIS_CARTESIAN)
      integer(c_int), value :: basis_type
      ! number of centers (atoms)
      integer(c_int), value :: num_centers
      ! center coordinates, dimension is 3*num_centers
      real(c_double)        :: center_xyz(*)
      ! center element (1 for H, 2 for He, ...), dimension is num_centers
      integer(c_int)        :: center_element(*)
      ! number of shells
      integer(c_int), value :: num_shells
      ! shell center, dimension is num_shells
      integer(c_int)        :: shell_center(*)
      ! L quantum number for each shell, dimension is num_shells
      integer(c_int)        :: l_quantum_num(*)
      ! number of primitives for each shell, dimension is num_shells
      integer(c_int)        :: shell_num_primitives(*)
      ! primitive exponents, dimension is sum(shell_num_primitives)
      real(c_double)        :: primitive_exp(*)
      ! contraction coefficients, dimension is sum(shell_num_primitives)
      real(c_double)        :: contraction_coef(*)
   end subroutine

   subroutine xcint_integrate(mode,          &
                              num_pert,      &
                              pert,          &
                              comp,          &
                              num_dmat,      &
                              dmat_to_pert,  &
                              dmat_to_comp,  &
                              dmat,          &
                              get_xc_energy, &
                              xc_energy,     &
                              get_xc_mat,    &
                              xc_mat,        &
                              num_electrons) &
              bind(c)

      integer(c_int), intent(in), value :: mode
      integer(c_int), intent(in), value :: num_pert
      integer(c_int), intent(in)        :: pert(*)
      integer(c_int), intent(in)        :: comp(*)
      integer(c_int), intent(in), value :: num_dmat
      integer(c_int), intent(in)        :: dmat_to_pert(*)
      integer(c_int), intent(in)        :: dmat_to_comp(*)
      real(c_double), intent(in)        :: dmat(*)
      integer(c_int), intent(in), value :: get_xc_energy
      real(c_double), intent(out)       :: xc_energy(*)
      integer(c_int), intent(in), value :: get_xc_mat
      real(c_double), intent(out)       :: xc_mat(*)
      real(c_double), intent(out)       :: num_electrons
   end subroutine

   subroutine xcint_set_verbosity(v) bind(c)
   ! input
      integer(c_int), intent(in), value :: v
   end subroutine

   subroutine xcint_set_stdout_function(fun) bind(c)
      type(c_funptr), intent(in), value :: fun
   ! input
      ! function with following signature
      ! integer(c_int) function fun(string) bind(c)
      !    character(kind=c_char, len=1), intent(in) :: string(*)
      ! end function
   end subroutine

   subroutine xcint_set_stderr_function(fun) bind(c)
      type(c_funptr), intent(in), value :: fun
   ! input
      ! function with following signature
      ! integer(c_int) function fun(string) bind(c)
      !    character(kind=c_char, len=1), intent(in) :: string(*)
      ! end function
   end subroutine

   subroutine xcint_integrate_worker() bind (c)
   end subroutine

   subroutine xcint_print_splash() bind (c)
   end subroutine

end module
