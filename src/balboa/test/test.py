import balboa
import numpy as np
from cffi import FFI
import os
import random


def sub(num_points,
        num_points_reference,
        generate_reference=False):

    assert num_points <= num_points_reference
    max_geo_order = 2
    num_slices = 10

    num_centers = 2
    center_coordinates_bohr = [
        1.7,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]

    num_shells = 9
    shell_centers = [1, 1, 1, 1, 1, 1, 2, 2, 2]
    shell_l_quantum_numbers = [0, 0, 0, 1, 1, 2, 0, 0, 1]
    shell_num_primitives = [9, 9, 1, 4, 1, 1, 4, 1, 1]

    primitive_exponents = [
        1.471e+04,
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
        7.270e-01,
    ]

    contraction_coefficients = [
        6.86365e-01,
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
        9.56881e-01,
    ]

    context = balboa.new_context()

    ierr = balboa.set_basis(context,
                            0,
                            num_centers,
                            center_coordinates_bohr,
                            num_shells,
                            shell_centers,
                            shell_l_quantum_numbers,
                            shell_num_primitives,
                            primitive_exponents,
                            contraction_coefficients)

    num_aos = balboa.get_num_aos(context)
    assert num_aos == 19

    dir_path = os.path.dirname(os.path.realpath(__file__))

    if generate_reference:
        random.seed(1)
        x_coordinates_bohr = []
        y_coordinates_bohr = []
        z_coordinates_bohr = []
        s = []
        for _ in range(num_points_reference):
            x = random.uniform(-2.0, 2.0)
            y = random.uniform(-2.0, 2.0)
            z = random.uniform(-2.0, 2.0)
            x_coordinates_bohr.append(x)
            y_coordinates_bohr.append(y)
            z_coordinates_bohr.append(z)
            s.append('{0} {1} {2}'.format(x, y, z))
        with open(os.path.join(dir_path, 'coordinates.txt'), 'w') as f:
            f.write('\n'.join(s))
    else:
        with open(os.path.join(dir_path, 'coordinates.txt'), 'r') as f:
            x_coordinates_bohr = []
            y_coordinates_bohr = []
            z_coordinates_bohr = []
            for line in f.readlines():
                (x, y, z) = line.split()
                x_coordinates_bohr.append(float(x))
                y_coordinates_bohr.append(float(y))
                z_coordinates_bohr.append(float(z))
        with open(os.path.join(dir_path, 'result.txt'), 'r') as f:
            ref_aos = []
            for line in f.readlines():
                ref_aos.append(float(line))

    # buffer length is adjusted for number of cartesian aos
    # so possibly longer than the number of spherical aos
    _l = balboa.get_buffer_len(context, max_geo_order, num_points)

    # allocate a numpy array of length l and zero it out
    aos = np.zeros(_l, dtype=np.float64)

    # we set the array to some number to make sure it is zeroed out inside
    aos.fill(123.456)

    # cast a pointer which points to the numpy array data
    ffi = FFI()
    aos_p = ffi.cast("double *", aos.ctypes.data)

    ierr = balboa.get_ao(context,
                         max_geo_order,
                         num_points,
                         x_coordinates_bohr,
                         y_coordinates_bohr,
                         z_coordinates_bohr,
                         aos_p)

    ao_centers = [balboa.get_ao_center(context, i) for i in range(num_aos)]
    assert ao_centers == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]

    assert balboa.get_geo_offset(context, 1, 0, 0) == 19
    assert balboa.get_geo_offset(context, 0, 1, 0) == 19 * 2
    assert balboa.get_geo_offset(context, 0, 0, 1) == 19 * 3
    assert balboa.get_geo_offset(context, 2, 0, 0) == 19 * 4
    assert balboa.get_geo_offset(context, 1, 1, 0) == 19 * 5
    assert balboa.get_geo_offset(context, 1, 0, 1) == 19 * 6
    assert balboa.get_geo_offset(context, 0, 2, 0) == 19 * 7
    assert balboa.get_geo_offset(context, 0, 1, 1) == 19 * 8

    if generate_reference:
        with open(os.path.join(dir_path, 'result.txt'), 'w') as f:
            k = 0
            for _diff in range(num_slices):
                for _ao in range(num_aos):
                    for _point in range(num_points):
                        f.write("{0}\n".format(aos_p[k]))
                        k += 1
    else:
        k = 0
        kr = 0
        for _diff in range(num_slices):
            for _ao in range(num_aos):
                for _point in range(num_points):
                    error = aos[k] - ref_aos[kr]
                    if abs(ref_aos[kr]) > 1.0e-20:
                        error /= ref_aos[kr]
                    assert abs(error) < 1.0e-8
                    kr += 1
                    k += 1
                kr += num_points_reference - num_points

    balboa.free_context(context)


def test_32():
    sub(num_points=32,
        num_points_reference=33)


def test_33():
    sub(num_points=33,
        num_points_reference=33)


if __name__ == '__main__':
    sub(num_points=33,
        num_points_reference=33,
        generate_reference=True)
