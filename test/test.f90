program test

   use xcint, only: xcint_new_context,    &
                    xcint_free_context,   &
                    xcint_set_basis,      &
                    xcint_set_functional, &
                    xcint_integrate_scf,  &
                    xcint_integrate,      &
                    XCINT_MODE_RKS,       &
                    XCINT_BASIS_SPHERICAL

   use numgrid

   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char

   implicit none

   type(c_ptr) :: numgrid_context
   type(c_ptr) :: xcint_context

   real(8)              :: radial_precision
   integer              :: min_num_angular_points
   integer              :: max_num_angular_points
   integer              :: num_centers
   real(8), allocatable :: center_coordinates(:)
   integer, allocatable :: proton_charges(:)
   integer              :: num_shells
   integer, allocatable :: shell_centers(:)
   integer, allocatable :: shell_l_quantum_numbers(:)
   integer, allocatable :: shell_num_primitives(:)
   real(8), allocatable :: primitive_exponents(:)
   real(8), allocatable :: contraction_coefficients(:)
   integer              :: num_points
   integer              :: num_points_center
   integer, parameter   :: io_unit = 13
   real(8)              :: ref(4)
   integer              :: i, j, k
   integer              :: center_index
   integer              :: ipoint
   real(8)              :: error
   integer              :: ierr
   real(8), allocatable :: dmat(:)
   real(8), allocatable :: vxc(:)
   integer              :: mat_dim
   real(8)              :: exc
   real(8)              :: num_electrons
   real(8)              :: f
   real(8)              :: alpha_max
   integer              :: max_l_quantum_number
   real(8)              :: alpha_min(3)
   real(8), allocatable :: atom_grid_x_bohr(:)
   real(8), allocatable :: atom_grid_y_bohr(:)
   real(8), allocatable :: atom_grid_z_bohr(:)
   real(8), allocatable :: atom_grid_w(:)
   real(8), allocatable :: grid_x_bohr(:)
   real(8), allocatable :: grid_y_bohr(:)
   real(8), allocatable :: grid_z_bohr(:)
   real(8), allocatable :: grid_w(:)
   real(8)              :: x_coordinates_bohr(2)
   real(8)              :: y_coordinates_bohr(2)
   real(8)              :: z_coordinates_bohr(2)

   radial_precision = 1.0d-12
   min_num_angular_points = 86
   max_num_angular_points = 302

   num_centers = 2

   allocate(center_coordinates(num_centers*3))

   center_coordinates(1) = 1.7d0
   center_coordinates(2) = 0.0d0
   center_coordinates(3) = 0.0d0
   center_coordinates(4) = 0.0d0
   center_coordinates(5) = 0.0d0
   center_coordinates(6) = 0.0d0

   x_coordinates_bohr(1) = 1.7d0
   x_coordinates_bohr(2) = 0.0d0

   y_coordinates_bohr(1) = 0.0d0
   y_coordinates_bohr(2) = 0.0d0

   z_coordinates_bohr(1) = 0.0d0
   z_coordinates_bohr(2) = 0.0d0

   allocate(proton_charges(num_centers))

   proton_charges(1) = 9
   proton_charges(2) = 1

   num_shells = 9

   allocate(shell_centers(num_shells))

   shell_centers(1) = 1
   shell_centers(2) = 1
   shell_centers(3) = 1
   shell_centers(4) = 1
   shell_centers(5) = 1
   shell_centers(6) = 1
   shell_centers(7) = 2
   shell_centers(8) = 2
   shell_centers(9) = 2

   allocate(shell_l_quantum_numbers(num_shells))

   shell_l_quantum_numbers(1) = 0
   shell_l_quantum_numbers(2) = 0
   shell_l_quantum_numbers(3) = 0
   shell_l_quantum_numbers(4) = 1
   shell_l_quantum_numbers(5) = 1
   shell_l_quantum_numbers(6) = 2
   shell_l_quantum_numbers(7) = 0
   shell_l_quantum_numbers(8) = 0
   shell_l_quantum_numbers(9) = 1

   allocate(shell_num_primitives(num_shells))

   shell_num_primitives(1) = 9
   shell_num_primitives(2) = 9
   shell_num_primitives(3) = 1
   shell_num_primitives(4) = 4
   shell_num_primitives(5) = 1
   shell_num_primitives(6) = 1
   shell_num_primitives(7) = 4
   shell_num_primitives(8) = 1
   shell_num_primitives(9) = 1

   allocate(primitive_exponents(31))

   primitive_exponents( 1) = 1.471d+4
   primitive_exponents( 2) = 2.207d+3
   primitive_exponents( 3) = 5.028d+2
   primitive_exponents( 4) = 1.426d+2
   primitive_exponents( 5) = 4.647d+1
   primitive_exponents( 6) = 1.670d+1
   primitive_exponents( 7) = 6.356d+0
   primitive_exponents( 8) = 1.316d+0
   primitive_exponents( 9) = 3.897d-1
   primitive_exponents(10) = 1.471d+4
   primitive_exponents(11) = 2.207d+3
   primitive_exponents(12) = 5.028d+2
   primitive_exponents(13) = 1.426d+2
   primitive_exponents(14) = 4.647d+1
   primitive_exponents(15) = 1.670d+1
   primitive_exponents(16) = 6.356d+0
   primitive_exponents(17) = 1.316d+0
   primitive_exponents(18) = 3.897d-1
   primitive_exponents(19) = 3.897d-1
   primitive_exponents(20) = 2.267d+1
   primitive_exponents(21) = 4.977d+0
   primitive_exponents(22) = 1.347d+0
   primitive_exponents(23) = 3.471d-1
   primitive_exponents(24) = 3.471d-1
   primitive_exponents(25) = 1.640d+0
   primitive_exponents(26) = 1.301d+1
   primitive_exponents(27) = 1.962d+0
   primitive_exponents(28) = 4.446d-1
   primitive_exponents(29) = 1.220d-1
   primitive_exponents(30) = 1.220d-1
   primitive_exponents(31) = 7.270d-1

   allocate(contraction_coefficients(31))

   contraction_coefficients( 1) =  6.86365d-1
   contraction_coefficients( 2) =  1.27435d+0
   contraction_coefficients( 3) =  2.13913d+0
   contraction_coefficients( 4) =  3.13055d+0
   contraction_coefficients( 5) =  3.63823d+0
   contraction_coefficients( 6) =  2.64148d+0
   contraction_coefficients( 7) =  7.55357d-1
   contraction_coefficients( 8) =  1.34270d-2
   contraction_coefficients( 9) = -8.19760d-4
   contraction_coefficients(10) = -1.57074d-1
   contraction_coefficients(11) = -3.00172d-1
   contraction_coefficients(12) = -4.91514d-1
   contraction_coefficients(13) = -7.84991d-1
   contraction_coefficients(14) = -9.34756d-1
   contraction_coefficients(15) = -1.00548d+0
   contraction_coefficients(16) = -3.20466d-1
   contraction_coefficients(17) =  4.92853d-1
   contraction_coefficients(18) =  1.99941d-1
   contraction_coefficients(19) =  3.51526d-1
   contraction_coefficients(20) =  3.16438d+0
   contraction_coefficients(21) =  2.49771d+0
   contraction_coefficients(22) =  1.05186d+0
   contraction_coefficients(23) =  1.73975d-1
   contraction_coefficients(24) =  3.79759d-1
   contraction_coefficients(25) =  6.77559d+0
   contraction_coefficients(26) =  9.61066d-2
   contraction_coefficients(27) =  1.63020d-1
   contraction_coefficients(28) =  1.85545d-1
   contraction_coefficients(29) =  7.37438d-2
   contraction_coefficients(30) =  1.47123d-1
   contraction_coefficients(31) =  9.56881d-1

   num_points = 31424  ! cheat to avoid reading twice
   allocate(grid_x_bohr(num_points))
   allocate(grid_y_bohr(num_points))
   allocate(grid_z_bohr(num_points))
   allocate(grid_w(num_points))

   num_points = 0
   do center_index = 1, num_centers

      ! ugly hack
      if (center_index == 1) then
         alpha_max = 14710.0d0
         max_l_quantum_number = 2
         alpha_min(1) = 0.3897d0
         alpha_min(2) = 0.3471d0
         alpha_min(3) = 1.64d0
      else
         alpha_max = 13.01d0
         max_l_quantum_number = 1
         alpha_min(1) = 0.122d0
         alpha_min(2) = 0.727d0
         alpha_min(3) = 0.0d0  ! not used
      end if

      numgrid_context = numgrid_new_atom_grid(radial_precision,             &
                                              min_num_angular_points,       &
                                              max_num_angular_points,       &
                                              proton_charges(center_index), &
                                              alpha_max,                    &
                                              max_l_quantum_number,         &
                                              alpha_min)

      num_points_center = numgrid_get_num_grid_points(numgrid_context)

      allocate(atom_grid_x_bohr(num_points_center))
      allocate(atom_grid_y_bohr(num_points_center))
      allocate(atom_grid_z_bohr(num_points_center))
      allocate(atom_grid_w(num_points_center))

      call numgrid_get_grid(numgrid_context,  &
                            num_centers,      &
                            center_index - 1, &
                            x_coordinates_bohr, &
                            y_coordinates_bohr, &
                            z_coordinates_bohr, &
                            proton_charges,   &
                            atom_grid_x_bohr,        &
                            atom_grid_y_bohr,        &
                            atom_grid_z_bohr,        &
                            atom_grid_w)

      do ipoint = 1, num_points_center
          grid_x_bohr(num_points + ipoint) = atom_grid_x_bohr(ipoint)
          grid_y_bohr(num_points + ipoint) = atom_grid_y_bohr(ipoint)
          grid_z_bohr(num_points + ipoint) = atom_grid_z_bohr(ipoint)
          grid_w(num_points + ipoint) = atom_grid_w(ipoint)
      end do
      num_points = num_points + num_points_center

      deallocate(atom_grid_x_bohr)
      deallocate(atom_grid_y_bohr)
      deallocate(atom_grid_z_bohr)
      deallocate(atom_grid_w)

      call numgrid_free_atom_grid(numgrid_context)
   end do

   xcint_context = xcint_new_context()

   ierr = xcint_set_basis(xcint_context,           &
                          XCINT_BASIS_SPHERICAL,   &
                          num_centers,             &
                          center_coordinates,      &
                          num_shells,              &
                          shell_centers,           &
                          shell_l_quantum_numbers, &
                          shell_num_primitives,    &
                          primitive_exponents,     &
                          contraction_coefficients)

   deallocate(center_coordinates)
   deallocate(proton_charges)
   deallocate(shell_centers)
   deallocate(shell_l_quantum_numbers)
   deallocate(shell_num_primitives)
   deallocate(primitive_exponents)
   deallocate(contraction_coefficients)

   mat_dim = 19

   allocate(dmat(mat_dim*mat_dim))
   dmat = 0.0d0
   open(unit=io_unit, file='../test/dmat.txt', access='sequential', action='read')
   do while (.true.)
      read(io_unit, *, end=10) k, f
      dmat(k+1) = f
   end do
10 close(io_unit)

   allocate(vxc(mat_dim*mat_dim))

   ierr = xcint_set_functional(xcint_context, "lda"//c_null_char)

   exc = 0.0
   num_electrons = 0.0

   ierr = xcint_integrate_scf(xcint_context,  &
                              XCINT_MODE_RKS, &
                              num_points,     &
                              grid_x_bohr,    &
                              grid_y_bohr,    &
                              grid_z_bohr,    &
                              grid_w,         &
                              dmat,           &
                              exc,            &
                              vxc,            &
                              num_electrons)

   if (dabs(num_electrons - 9.999992072209077d0) > 1.0e-12) stop 1
   if (dabs(exc + 20.421064966253642d0) > 1.0e-12) stop 1

   ierr = xcint_set_functional(xcint_context, "b3lyp"//c_null_char)

   ierr = xcint_integrate(xcint_context,  &
                          XCINT_MODE_RKS, &
                          num_points,     &
                          grid_x_bohr,    &
                          grid_y_bohr,    &
                          grid_z_bohr,    &
                          grid_w,         &
                          0,              &
                          (/0/),          &
                          (/0/),          &
                          0,              &
                          (/0/),          &
                          dmat,           &
                          1,              &
                          exc,            &
                          1,              &
                          vxc,            &
                          num_electrons)

   deallocate(grid_x_bohr)
   deallocate(grid_y_bohr)
   deallocate(grid_z_bohr)
   deallocate(grid_w)

   if (dabs(num_electrons - 9.999992072209077d0) > 1.0e-12) stop 1
   if (dabs(exc + 17.475254754225027d0) > 1.0e-12) stop 1

   deallocate(dmat)
   deallocate(vxc)

   call xcint_free_context(xcint_context)

end program
