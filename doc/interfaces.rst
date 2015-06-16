

==========
Interfaces
==========


General
=======

Functions that interact with user return 0 upon success.
Programming errors fail immediately and do not return.


xcint_set_functional
====================

Set the XC functional.

.. code-block:: fortran

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


xcint_generate_grid
===================

Generate the integration grid.

.. code-block:: fortran

   integer(c_int) function xcint_generate_grid(radial_precision,     &
                                               angular_min,          &
                                               angular_max,          &
                                               num_centers,          &
                                               center_coordinates,           &
                                               center_elements,       &
                                               num_shells,           &
                                               shell_centers,         &
                                               l_quantum_numbers,        &
                                               shell_num_primitives, &
                                               primitive_exponents)        &
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
      real(c_double)        :: center_coordinates(*)
      ! center element (1 for H, 2 for He, ...), dimension is num_centers
      integer(c_int)        :: center_elements(*)
      ! number of shells
      integer(c_int), value :: num_shells
      ! shell center, dimension is num_shells
      integer(c_int)        :: shell_centers(*)
      ! L quantum number for each shell, dimension is num_shells
      integer(c_int)        :: l_quantum_numbers(*)
      ! number of primitives for each shell, dimension is num_shells
      integer(c_int)        :: shell_num_primitives(*)
      ! primitive exponents, dimension is sum(shell_num_primitives)
      real(c_double)        :: primitive_exponents(*)
   ! returns
   !    0 upon success
   end function


xcint_set_basis
===============

Set the basis set.

.. code-block:: fortran

   subroutine xcint_set_basis(basis_type,           &
                              num_centers,          &
                              center_coordinates,           &
                              center_elements,       &
                              num_shells,           &
                              shell_centers,         &
                              l_quantum_numbers,        &
                              shell_num_primitives, &
                              primitive_exponents,        &
                              contraction_coefficients)     &
              bind (c)
   ! input
      ! basis set type (XCINT_BASIS_SPHERICAL or XCINT_BASIS_CARTESIAN)
      integer(c_int), value :: basis_type
      ! number of centers (atoms)
      integer(c_int), value :: num_centers
      ! center coordinates, dimension is 3*num_centers
      real(c_double)        :: center_coordinates(*)
      ! center element (1 for H, 2 for He, ...), dimension is num_centers
      integer(c_int)        :: center_elements(*)
      ! number of shells
      integer(c_int), value :: num_shells
      ! shell center, dimension is num_shells
      integer(c_int)        :: shell_centers(*)
      ! L quantum number for each shell, dimension is num_shells
      integer(c_int)        :: l_quantum_numbers(*)
      ! number of primitives for each shell, dimension is num_shells
      integer(c_int)        :: shell_num_primitives(*)
      ! primitive exponents, dimension is sum(shell_num_primitives)
      real(c_double)        :: primitive_exponents(*)
      ! contraction coefficients, dimension is sum(shell_num_primitives)
      real(c_double)        :: contraction_coefficients(*)
   end subroutine


xcint_integrate
===============

The workhorse of XCint: integrate XC energies and matrix elements.

.. code-block:: fortran

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


Arguments
---------

**mode** (input)

Possible entries:

- XCINT_MODE_RKS -- Restricted Kohn-Sham.
- XCINT_MODE_UKS -- Unrestricted Kohn-Sham (currently not supported).


**num_pert** (input)

Number of perturbations. Has to be 0 or positive integer.


**pert** (input)

Dimension is num_pert.

Not used if num_pert is 0.

Possible perturbation types:

- XCINT_PERT_EL -- Electric perturbation.
- XCINT_PERT_GEO -- Geometric perturbation.
- XCINT_PERT_MAG_CGO -- Magnetic perturbation (currently not supported).
- XCINT_PERT_MAG_LAO -- London AO magnetic perturbation (currently not supported).


**comp** (input)

Dimension is 2*num_pert.

Not used if num_pert is 0.

For each perturbation we expect 2 integers:
start component and end component.


**num_dmat** (input)

Number of density matrices.


**dmat_to_pert** (input)

Mapping of density matrices to perturbations.


**dmat_to_comp** (input)

Mapping of density matrices to perturbation components.


**dmat** (input)

Array that holds the density matrix or matrices.


**get_xc_energy** (input)

- 0 -- Do not integrate the XC energy (derivatives).
- 1 -- Integrate the XC energy (derivatives).


**xc_energy** (output)

Array that holds the integrated XC energy (derivative or derivatives).

Not touched if get_xc_energy is 0.


**get_xc_mat** (input)

- 0 -- Do not integrate the XC potential matrix (derivatives).
- 1 -- Integrate the XC potential matrix (derivatives).


**xc_mat** (output)

Array that holds the integrated XC potential matrix (derivative or derivatives).

Not touched if get_xc_mat is 0.


**num_electrons** (output)

Integrated number of electrons.


xcint_set_stdout_function
=========================

Function which implements printing to "stdout".

.. code-block:: fortran

   subroutine xcint_set_stdout_function(fun) bind(c)
      type(c_funptr), intent(in), value :: fun
   ! input
      ! function with following signature
      ! integer(c_int) function fun(string) bind(c)
      !    character(kind=c_char, len=1), intent(in) :: string(*)
      ! end function
   end subroutine


xcint_set_stderr_function
=========================

Function which implements printing to "stderr".

.. code-block:: fortran

   subroutine xcint_set_stderr_function(fun) bind(c)
      type(c_funptr), intent(in), value :: fun
   ! input
      ! function with following signature
      ! integer(c_int) function fun(string) bind(c)
      !    character(kind=c_char, len=1), intent(in) :: string(*)
      ! end function
   end subroutine


xcint_integrate_worker
======================

Starts the MPI worker process.

.. code-block:: fortran

   subroutine xcint_integrate_worker() bind (c)
   end subroutine


xcint_print_splash
==================

Print splash screen

.. code-block:: fortran

   subroutine xcint_print_splash() bind (c)
   end subroutine
