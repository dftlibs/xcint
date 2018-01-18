module xcint

   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int, c_char

   implicit none

   public xcint_new_context
   public xcint_free_context
   public xcint_set_functional
   public xcint_set_basis
   public xcint_integrate_scf
   public xcint_integrate

   public XCINT_MODE_RKS
   public XCINT_MODE_UKS
   public XCINT_BASIS_SPHERICAL
   public XCINT_BASIS_CARTESIAN
   public XCINT_PERT_EL
   public XCINT_PERT_GEO
   public XCINT_PERT_MAG_CGO
   public XCINT_PERT_MAG_LAO

   private

   enum, bind(c)
      enumerator :: XCINT_MODE_RKS
      enumerator :: XCINT_MODE_UKS
   end enum

   enum, bind(c)
      enumerator :: XCINT_BASIS_SPHERICAL
      enumerator :: XCINT_BASIS_CARTESIAN
   end enum

   enum, bind(c)
      enumerator :: XCINT_PERT_EL
      enumerator :: XCINT_PERT_GEO
      enumerator :: XCINT_PERT_MAG_CGO
      enumerator :: XCINT_PERT_MAG_LAO
   end enum

   interface xcint_new_context
      function xcint_new_context() result(context) bind (C)
         import :: c_ptr
         type(c_ptr) :: context
      end function
   end interface

   interface xcint_free_context
      subroutine xcint_free_context(context) bind (C)
         import :: c_ptr
         type(c_ptr), value :: context
      end subroutine
   end interface

   interface xcint_set_functional
      function xcint_set_functional(context,                 &
                                    line) result(ierr) bind (C)
         import :: c_ptr, c_int, c_char
         type(c_ptr), value            :: context
         character(c_char), intent(in) :: line
         integer(c_int) :: ierr
      end function
   end interface

   interface xcint_set_basis
      function xcint_set_basis(context,                 &
                               basis_type,              &
                               num_centers,             &
                               center_coordinates,      &
                               num_shells,              &
                               shell_centers,           &
                               shell_l_quantum_numbers, &
                               shell_num_primitives,    &
                               primitive_exponents,     &
                               contraction_coefficients) result(ierr) bind (C)
         import :: c_ptr, c_int, c_double
         type(c_ptr), value                :: context
         integer(c_int), intent(in), value :: basis_type
         integer(c_int), intent(in), value :: num_centers
         real(c_double), intent(in)        :: center_coordinates(*)
         integer(c_int), intent(in), value :: num_shells
         integer(c_int), intent(in)        :: shell_centers(*)
         integer(c_int), intent(in)        :: shell_l_quantum_numbers(*)
         integer(c_int), intent(in)        :: shell_num_primitives(*)
         real(c_double), intent(in)        :: primitive_exponents(*)
         real(c_double), intent(in)        :: contraction_coefficients(*)
         integer(c_int) :: ierr
      end function
   end interface

   interface xcint_integrate_scf
      function xcint_integrate_scf(context,       &
                                   mode,          &
                                   num_points,    &
                                   grid_x_bohr,   &
                                   grid_y_bohr,   &
                                   grid_z_bohr,   &
                                   grid_w,        &
                                   dmat,          &
                                   exc,           &
                                   vxc,           &
                                   num_electrons) result(ierr) bind (C)
         import :: c_ptr, c_int, c_double
         type(c_ptr), value                :: context
         integer(c_int), intent(in), value :: mode
         integer(c_int), intent(in), value :: num_points
         real(c_double), intent(in)        :: grid_x_bohr(*)
         real(c_double), intent(in)        :: grid_y_bohr(*)
         real(c_double), intent(in)        :: grid_z_bohr(*)
         real(c_double), intent(in)        :: grid_w(*)
         real(c_double), intent(in)        :: dmat(*)
         real(c_double), intent(inout)     :: exc
         real(c_double), intent(inout)     :: vxc(*)
         real(c_double), intent(inout)     :: num_electrons
         integer(c_int) :: ierr
      end function
   end interface

   interface xcint_integrate
      function xcint_integrate(context,              &
                               mode,                 &
                               num_points,           &
                               grid_x_bohr,          &
                               grid_y_bohr,          &
                               grid_z_bohr,          &
                               grid_w,               &
                               num_perturbations,    &
                               perturbations,        &
                               components,           &
                               num_dmat,             &
                               perturbation_indices, &
                               dmat,                 &
                               get_exc,              &
                               exc,                  &
                               get_vxc,              &
                               vxc,                  &
                               num_electrons) result(ierr) bind (C)
         import :: c_ptr, c_int, c_double
         type(c_ptr), value                :: context
         integer(c_int), intent(in), value :: mode
         integer(c_int), intent(in), value :: num_points
         real(c_double), intent(in)        :: grid_x_bohr(*)
         real(c_double), intent(in)        :: grid_y_bohr(*)
         real(c_double), intent(in)        :: grid_z_bohr(*)
         real(c_double), intent(in)        :: grid_w(*)
         integer(c_int), intent(in), value :: num_perturbations
         integer(c_int), intent(in)        :: perturbations(*)
         integer(c_int), intent(in)        :: components(*)
         integer(c_int), intent(in), value :: num_dmat
         integer(c_int), intent(in)        :: perturbation_indices(*)
         real(c_double), intent(in)        :: dmat(*)
         integer(c_int), intent(in), value :: get_exc
         real(c_double), intent(inout)     :: exc
         integer(c_int), intent(in), value :: get_vxc
         real(c_double), intent(inout)     :: vxc(*)
         real(c_double), intent(inout)     :: num_electrons
         integer(c_int) :: ierr
      end function
   end interface

end module
