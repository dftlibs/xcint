

import sys


def get_c_code(l):
    s = ''
    for i, x in enumerate(l):
        s += "const int %s = %i;\n" % (x, i)
    return s


def get_fortran_code(l):
    s = ''
    for i, x in enumerate(l):
        s += "integer(c_int), parameter :: %s = %i\n" % (x, i)
    return s


def main(argv):

    l = []

    l.append('XCINT_BASIS_SPHERICAL')
    l.append('XCINT_BASIS_CARTESIAN')
    l.append('XCINT_MODE_RKS')
    l.append('XCINT_MODE_UKS')
    l.append('XCINT_PERT_EL')
    l.append('XCINT_PERT_GEO')
    l.append('XCINT_PERT_MAG_CGO')
    l.append('XCINT_PERT_MAG_LAO')

    f = open(argv[1], 'w')
    f.write(get_c_code(l))
    f.close()

    f = open(argv[2], 'w')
    f.write(get_fortran_code(l))
    f.close()

if __name__ == '__main__':
    main(sys.argv)
