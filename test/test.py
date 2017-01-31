import pytest


def test_energy():
    """
    Test energy.
    """
    import os
    import xcint
    import numgrid

    radial_precision = 1.0e-12
    min_num_angular_points = 86
    max_num_angular_points = 302

    num_centers = 2
    center_coordinates = [1.7,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0]
    center_elements = [9, 1]

    num_outer_centers = 0
    outer_center_coordinates = []
    outer_center_elements = []

    num_shells = 9
    shell_centers = [1, 1, 1, 1, 1, 1, 2, 2, 2]
    shell_l_quantum_numbers = [0, 0, 0, 1, 1, 2, 0, 0, 1]
    shell_num_primitives = [9, 9, 1, 4, 1, 1, 4, 1, 1]

    primitive_exponents = [1.471e+04,
                           2.207e+03,
                           5.028e+02,
                           1.426e+02,
                           4.647e+01,
                           1.670e+01,
                           6.356e+00,
                           1.316e+00,
                           3.897e-01,
                           1.471e+04,
                           2.207e+03,
                           5.028e+02,
                           1.426e+02,
                           4.647e+01,
                           1.670e+01,
                           6.356e+00,
                           1.316e+00,
                           3.897e-01,
                           3.897e-01,
                           2.267e+01,
                           4.977e+00,
                           1.347e+00,
                           3.471e-01,
                           3.471e-01,
                           1.640e+00,
                           1.301e+01,
                           1.962e+00,
                           4.446e-01,
                           1.220e-01,
                           1.220e-01,
                           7.270e-01]

    numgrid_context = numgrid.new_context()

    ierr = numgrid.generate_grid(numgrid_context,
                                 radial_precision,
                                 min_num_angular_points,
                                 max_num_angular_points,
                                 num_centers,
                                 center_coordinates,
                                 center_elements,
                                 num_outer_centers,
                                 outer_center_coordinates,
                                 outer_center_elements,
                                 num_shells,
                                 shell_centers,
                                 shell_l_quantum_numbers,
                                 shell_num_primitives,
                                 primitive_exponents)

    num_points = numgrid.get_num_points(numgrid_context)

    assert num_points == 31424

    grid = numgrid.get_grid(numgrid_context)

    numgrid.free_context(numgrid_context)

    xcint_context = xcint.lib.xcint_new_context()

    contraction_coefficients = [6.86365e-01,
                                1.27435e+00,
                                2.13913e+00,
                                3.13055e+00,
                                3.63823e+00,
                                2.64148e+00,
                                7.55357e-01,
                                1.34270e-02,
                                -8.19760e-04,
                                -1.57074e-01,
                                -3.00172e-01,
                                -4.91514e-01,
                                -7.84991e-01,
                                -9.34756e-01,
                                -1.00548e+00,
                                -3.20466e-01,
                                4.92853e-01,
                                1.99941e-01,
                                3.51526e-01,
                                3.16438e+00,
                                2.49771e+00,
                                1.05186e+00,
                                1.73975e-01,
                                3.79759e-01,
                                6.77559e+00,
                                9.61066e-02,
                                1.63020e-01,
                                1.85545e-01,
                                7.37438e-02,
                                1.47123e-01,
                                9.56881e-01]

    ierr = xcint.lib.xcint_set_basis(xcint_context,
                                     xcint.lib.XCINT_BASIS_SPHERICAL,
                                     num_centers,
                                     center_coordinates,
                                     num_shells,
                                     shell_centers,
                                     shell_l_quantum_numbers,
                                     shell_num_primitives,
                                     primitive_exponents,
                                     contraction_coefficients)

    mat_dim = 19
    dmat = []
    for i in range(mat_dim * mat_dim):
        dmat.append(0.0)

    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, 'dmat.txt'), 'r') as f:
        for line in f.readlines():
            (i, d) = line.split()
            dmat[int(i)] = float(d)

    ierr = xcint.lib.xcint_set_functional(xcint_context, "lda")

    exc = xcint.ffi.new("double *")
    num_electrons = xcint.ffi.new("double *")
    vxc = xcint.ffi.new("double[]", mat_dim * mat_dim)

    ierr = xcint.lib.xcint_integrate_scf(xcint_context,
                                         xcint.lib.XCINT_MODE_RKS,
                                         num_points,
                                         grid,
                                         dmat,
                                         exc,
                                         vxc,
                                         num_electrons)

    assert abs(num_electrons[0] - 9.999992074832) < 1.0e-11
    assert abs(exc[0] - -20.421064966255539) < 1.0e-11

    dot = 0.0
    for i in range(mat_dim * mat_dim):
        dot += vxc[i] * dmat[i]
    assert abs(dot - -6.729996811122) < 1.0e-11

    xcint.lib.xcint_free_context(xcint_context)
