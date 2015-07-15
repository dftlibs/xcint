
import sys
import os
import cs_trans

AO_BLOCK_LENGTH = 128
AO_CHUNK_LENGTH = 32
MAX_GEO_DIFF_ORDER = 5
MAX_L_VALUE = 5

#-------------------------------------------------------------------------------

cs = cs_trans.get_cs_trans(MAX_L_VALUE)

#-------------------------------------------------------------------------------

def get_ijk_list(m):
    l = []
    for a in range(1, m + 2):
        for b in range(1, a + 1):
            i = m + 1 - a
            j = a - b
            k = b - 1
            l.append([i, j, k])
    return l

#-------------------------------------------------------------------------------

exp_offset = {}
for l in range(0, MAX_L_VALUE+1):
    k = 0
    for exp in get_ijk_list(l):
        exp_offset[tuple(exp)] = k
        k += 1

#-------------------------------------------------------------------------------

def get_exp_offset(exp):
    return exp_offset[tuple(exp)]

#-------------------------------------------------------------------------------

def get_ao_pointer_prefix(geo):
    return 'ao_%i%i%i' % (geo[0], geo[1], geo[2])

#-------------------------------------------------------------------------------

def write_parameters(file_name):
    f = open(file_name, 'w')
    f.write('const int AO_BLOCK_LENGTH = %i;\n' % AO_BLOCK_LENGTH)
    f.write('const int AO_CHUNK_LENGTH = %i;\n' % AO_CHUNK_LENGTH)
    f.write('const int MAX_GEO_DIFF_ORDER = %i;\n' % MAX_GEO_DIFF_ORDER)
    f.write('const int MAX_L_VALUE = %i;\n' % MAX_L_VALUE)
    f.close()

#-------------------------------------------------------------------------------

def print_line(exp, geo, m, r):

    exp_right = exp[:]
    exp_right[m] -= 1

    vec_r = 'buffer[OFFSET_%02d_%02d_%02d_%i%i%i]' % (exp[0], exp[1], exp[2], geo[0], geo[1], geo[2])
    vec_p = '%s' % r
    vec_a = 'buffer[OFFSET_%02d_%02d_%02d_%i%i%i]' % (exp_right[0], exp_right[1], exp_right[2], geo[0], geo[1], geo[2])

    if (geo[m] > 0):
        geo_right = geo[:]
        geo_right[m] -= 1
        vec_b = 'buffer[OFFSET_%02d_%02d_%02d_%i%i%i]' % (exp_right[0], exp_right[1], exp_right[2], geo_right[0], geo_right[1], geo_right[2])
        if geo[m] > 1:
            return '            get_pa_plus_sb(&%s, %s, %i.0, &%s, &%s);\n' % (vec_a, vec_p, geo[m], vec_b, vec_r)
        else:
            return '            get_pa_plus_b(&%s, %s, &%s, &%s);\n' % (vec_a, vec_p, vec_b, vec_r)
    else:
        return '            get_pa(&%s, %s, &%s);\n' % (vec_a, vec_p, vec_r)

#-------------------------------------------------------------------------------

def write_offsets(file_name):

    s = '#ifndef offsets_h_\n#define offsets_h_\n\n'
    offset = 0
    for l in range(0, MAX_L_VALUE+1):
        for exp in get_ijk_list(l):
            for g in range(0, MAX_GEO_DIFF_ORDER+1):
                for geo in get_ijk_list(g):
                    s_geo = '%i%i%i' % (geo[0], geo[1], geo[2])
                    s += '#define OFFSET_%02d_%02d_%02d_%s %i\n' % (exp[0], exp[1], exp[2], s_geo, offset)
                    offset += AO_CHUNK_LENGTH
    s += '\n#define BUFFER_LENGTH %i\n' % offset
    s += '\n#endif // offsets_h_\n'
    f = open(file_name, 'w')
    f.write(s)
    f.close()

#-------------------------------------------------------------------------------

def write_routine(_maxg, file_name):

    sfoo = '''
//  this file is automatically generated by generate.py

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <cstring>

#include "autogenerated.h"
#include "offsets.h"
#include "ao_vector.h"
          \n'''

    s = '''void get_ao_g%i(const int    shell_l_quantum_numbers,
               const int    num_primitives,
               const bool   is_spherical,
               const double primitive_exponents[],
               const double contraction_coefficients[],
                     double s[],
                     double buffer[],
               const double shell_centers_coordinates[],
               const double extent_squared,
               const double pw[],
                     double px[],
                     double py[],
                     double pz[],
                     double p2[],\n''' % _maxg

    l = []
    for g in range(0, _maxg+1):
        for geo in get_ijk_list(g):
            l.append('double ao_%i%i%i[]' % (geo[0], geo[1], geo[2]))

    for i in range(len(l)):
        if i < len(l)-1:
            s += '                     %s,\n' % l[i]
        else:
            s += '                     %s)' % l[i]

    sfoo += '''
%s
{
         \n''' % (s)

    if _maxg > 0:
        for i in range(_maxg + 1):
            sfoo += '    double fx_%i;\n' % i
            sfoo += '    double fy_%i;\n' % i
            sfoo += '    double fz_%i;\n' % i
    sfoo += '    double a;\n'
    sfoo += '    double c;\n\n'

    sfoo += '    for (int koff = 0; koff < %i; koff += %i)\n' % (AO_BLOCK_LENGTH, AO_CHUNK_LENGTH)
    sfoo += '    {\n'

    sfoo += '        get_p2(shell_centers_coordinates,\n'
    sfoo += '               &pw[4*koff],\n'
    sfoo += '               px,\n'
    sfoo += '               py,\n'
    sfoo += '               pz,\n'
    sfoo += '               p2);\n\n'

    sfoo += '        // screening\n'
    sfoo += '        if (not calculate_chunk(extent_squared, p2)) continue;\n\n'

    array = 'buffer[OFFSET_00_00_00_000]'
    sfoo += '        memset(&%s, 0, %i*sizeof(double));\n' % (array, AO_CHUNK_LENGTH)
    for g in range(1, _maxg+1):
        for geo in get_ijk_list(g):
            array = 'buffer[OFFSET_00_00_00_%i%i%i]' % (geo[0], geo[1], geo[2])
            sfoo += '        memset(&%s, 0, %i*sizeof(double));\n' % (array, AO_CHUNK_LENGTH)

    sfoo += '''
        for (int i = 0; i < num_primitives; i++)
        {
            a = -primitive_exponents[i];
            c = contraction_coefficients[i];

            get_exp(p2, c, a, s);

            #pragma ivdep
            #pragma vector aligned
            for (int k = 0; k < %i; k++)
            {
                buffer[OFFSET_00_00_00_000 + k] += s[k];
              \n''' % AO_CHUNK_LENGTH

    if (_maxg > 0):
        sfoo += '''                fx_0 = 1.0;
                fy_0 = 1.0;
                fz_0 = 1.0;
                fx_1 = 2.0*a*px[k];
                fy_1 = 2.0*a*py[k];
                fz_1 = 2.0*a*pz[k];
                buffer[OFFSET_00_00_00_100 + k] += fx_1*s[k];
                buffer[OFFSET_00_00_00_010 + k] += fy_1*s[k];
                buffer[OFFSET_00_00_00_001 + k] += fz_1*s[k];
                  \n'''

    for g in range(2, _maxg+1):
        sfoo += '                fx_%i = fx_%i*fx_1 + %i.0*a*fx_%i;\n' % (int(g), int(g-1), int(g-1)*2, int(g-2))
        sfoo += '                fy_%i = fy_%i*fy_1 + %i.0*a*fy_%i;\n' % (int(g), int(g-1), int(g-1)*2, int(g-2))
        sfoo += '                fz_%i = fz_%i*fz_1 + %i.0*a*fz_%i;\n' % (int(g), int(g-1), int(g-1)*2, int(g-2))
        for geo in get_ijk_list(g):
            sfoo += '                buffer[OFFSET_00_00_00_%i%i%i + k] += fx_%i*fy_%i*fz_%i*s[k];\n' \
                                                    % (geo[0], geo[1], geo[2], \
                                                       geo[0], geo[1], geo[2])
    sfoo += '            }\n'
    sfoo += '        }\n'

    for l in range(0, MAX_L_VALUE+1):
        sfoo += '\n        if (shell_l_quantum_numbers == ' + '%i)\n' % l
        sfoo += '        {\n'
        if l < 2:
            c = 0
            for exp in get_ijk_list(l):
                for s in range(len(cs[l][c])):
                    f = cs[l][c][s]
                    if abs(f) > 0.0:
                        for g in range(0, _maxg+1):
                            for geo in get_ijk_list(g):
                                s_geo = '%i%i%i' % (geo[0], geo[1], geo[2])
                                sfoo += '            memcpy(&%s[%i*%i + koff], &buffer[OFFSET_%02d_%02d_%02d_%s], %i*sizeof(double));\n' % (get_ao_pointer_prefix(geo), s, AO_BLOCK_LENGTH, exp[0], exp[1], exp[2], s_geo, AO_CHUNK_LENGTH)
                c += 1
        else:
            sfoo += '            if (is_spherical)\n'
            sfoo += '            {\n'
            c = 0
            for exp in get_ijk_list(l):
                for s in range(len(cs[l][c])):
                    f = cs[l][c][s]
                    if abs(f) > 0.0:
                        for g in range(0, _maxg+1):
                            for geo in get_ijk_list(g):
                                s_geo = '%i%i%i' % (geo[0], geo[1], geo[2])
                                sfoo += '                vec_daxpy(%20.16e, &buffer[OFFSET_%02d_%02d_%02d_%s], &%s[%i*%i + koff]);\n' % (f,
                                                                                                                                 exp[0], exp[1], exp[2], s_geo,
                                                                                                                                 get_ao_pointer_prefix(geo),
                                                                                                                                 s, AO_BLOCK_LENGTH)
                c += 1
            sfoo += '            }\n'
            sfoo += '            else\n'
            sfoo += '            {\n'
            s = 0
            for exp in get_ijk_list(l):
                for g in range(0, _maxg+1):
                    for geo in get_ijk_list(g):
                        s_geo = '%i%i%i' % (geo[0], geo[1], geo[2])
                        sfoo += '                memcpy(&%s[%i*%i + koff], &buffer[OFFSET_%02d_%02d_%02d_%s], %i*sizeof(double));\n' % (get_ao_pointer_prefix(geo), s, AO_BLOCK_LENGTH, exp[0], exp[1], exp[2], s_geo, AO_CHUNK_LENGTH)
                s += 1
            sfoo += '            }\n'
        sfoo += '            continue;\n'
        sfoo += '        }\n'
        sfoo += '        else\n'
        sfoo += '        {\n'
        if l+1 < MAX_L_VALUE+1:
            for exp in get_ijk_list(l+1):
                for g in range(0, _maxg+1):
                    for geo in get_ijk_list(g):
                        if exp[0] > 0:
                            sfoo += print_line(exp, geo, 0, 'px')
                        else:
                            if exp[1] > 0:
                                sfoo += print_line(exp, geo, 1, 'py')
                            else:
                                if exp[2] > 0:
                                    sfoo += print_line(exp, geo, 2, 'pz')
        else:
            sfoo += '             std::cout << "error: order too high";\n'
            sfoo += '             exit(1);\n'
        sfoo += '        }\n'
    sfoo += '    }\n'
    sfoo += '\n}\n'

    f = open(file_name, 'w')
    f.write(sfoo)
    f.close()

#-------------------------------------------------------------------------------

def write_aocalls(file_name):

    f = open(file_name, 'w')

    s1 = '''              basis.shell_l_quantum_numbers[ishell],
              basis.shell_num_primitives[ishell],
              basis.is_spherical,
              &basis.primitive_exponents[n],
              &basis.contraction_coefficients[n],
              s,
              buffer,
              &basis.shell_centers_coordinates[3*ishell],
              basis.shell_extent_squared[ishell],
              p,
              px,
              py,
              pz,
              p2,
              &ao_local[basis.shell_off[ishell]*%i]''' % AO_BLOCK_LENGTH
    s3 = '''
             );
    break;'''

    for g in range(0, MAX_GEO_DIFF_ORDER+1):

        f.write('\ncase %i:\n' % g)
        f.write('    get_ao_g%i(\n' % g)

        s2 = s1
        if g > 0:
            j = 0
            for _g in range(1, g+1):
                for geo in get_ijk_list(_g):
                    j += 1
                    s2 += ',\n              &ao_local[(basis.shell_off[ishell] + %i*basis.num_ao)*%i]' % (j, AO_BLOCK_LENGTH)

        f.write(s2)
        f.write(s3)

    f.close()

#-------------------------------------------------------------------------------

def write_header(file_name):

    f = open(file_name, 'w')

    f.write('#ifndef autogenerated_h_\n')
    f.write('#define autogenerated_h_\n')
    f.write('\nextern "C"\n')
    f.write('{\n')

    for g in range(0, MAX_GEO_DIFF_ORDER+1):
        f.write('    void get_ao_g%i(const int    shell_l_quantum_numbers,\n' % g)
        f.write('                   const int    num_primitives,\n')
        f.write('                   const bool   is_spherical,\n')
        f.write('                   const double primitive_exponents[],\n')
        f.write('                   const double contraction_coefficients[],\n')
        f.write('                         double s[],\n')
        f.write('                         double buffer[],\n')
        f.write('                   const double shell_centers_coordinates[],\n')
        f.write('                   const double extent_squared,\n')
        f.write('                   const double pw[],\n')
        f.write('                         double px[],\n')
        f.write('                         double py[],\n')
        f.write('                         double pz[],\n')
        f.write('                         double p2[]')

        s = ''
        for _g in range(g+1):
            for geo in get_ijk_list(_g):
                s += ',\n                         double ao_%i%i%i[]' % (geo[0], geo[1], geo[2])

        f.write(s)
        f.write('\n                  );\n')

    f.write('}\n\n')
    f.write('#endif // autogenerated_h_\n')
    f.close()

#-------------------------------------------------------------------------------

def main(output_directory):
    write_parameters(os.path.join(output_directory, 'parameters.h'))
    write_offsets(os.path.join(output_directory, 'offsets.h'))
    for g in range(0, MAX_GEO_DIFF_ORDER+1):
        write_routine(g, os.path.join(output_directory, 'autogenerated_%i.cpp' % g))
    write_aocalls(os.path.join(output_directory, 'aocalls.h'))
    write_header(os.path.join(output_directory, 'autogenerated.h'))

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv[1])
