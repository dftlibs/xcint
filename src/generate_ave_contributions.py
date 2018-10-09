MAX_ORDER = 6

import math
import itertools


def add_perturbation(l):
    order = int(math.log(len(l)+1, 2))
    l_new = []
    for x in l:
        l_new.append(x)
    l_new.append([order+1])
    for x in l:
        y = x[:]
        y.extend([order+1])
        l_new.append(y)
    return l_new


def get_fortran_code(j, density_index, dmat_index_rest):
    s  = ''
    s += 'if (geo_derv_order > %i)\n' % (max(j) - 1)
    s += '{\n'
    s += '    if (use_dmat[%i])\n' % dmat_index_rest
    s += '    {\n'
    s += '        if (!n_is_used[%i])\n' % density_index
    s += '        {\n'
    s += '            std::fill(&n[%i*block_length*num_variables], &n[(%i+1)*block_length*num_variables], 0.0);\n' % (density_index, density_index)
    s += '            n_is_used[%i] = true;\n' % density_index
    s += '        }\n'
    for i in j:
        s += '        coor.push_back(geo_coor[%i]);\n' % int(i-1)
    s += '            get_dens_geo_derv(mat_dim,\n'
    s += '                              num_aos,\n'
    s += '                              block_length,\n'
    s += '                              buffer_len,\n'
    s += '                              ao,\n'
    s += '                              ao_centers,\n'
    s += '                              get_gradient,\n'
    s += '                              get_tau,\n'
    s += '                              coor,\n'
    s += '                              get_geo_offset,\n'
    s += '                              &n[%i*block_length*num_variables],\n' % density_index
    s += '                              &dmat[dmat_index[%i]]);\n' % dmat_index_rest
    s += '        coor.clear();\n'
    s += '    }\n'
    s += '}\n'
    return s


def get_all_possible_selections(l, d):
    for i in range(len(l)):
        for j in itertools.combinations(l, i+1):
            k = l[:]
            for m in j:
                k.remove(m)
            density_index   = d[tuple(l)]
            dmat_index_rest = d[tuple(k)]
            print(get_fortran_code(j, density_index, dmat_index_rest))


def main():

    print('// this file is autogenerated by generate_ave_contributions.py\n')

    l = [[1]]
    for i in range(MAX_ORDER - 1):
        l = add_perturbation(l)

    d = {}
    d[()] = 0
    for i, x in enumerate(l):
        d[tuple(x)] = i+1

    for x in l:
        get_all_possible_selections(x, d)


if __name__ == '__main__':
    main()
