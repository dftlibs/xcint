
"""
    CMake-boot: bootstrapping front-end to CMake
    Copyright 2014 by Radovan Bast and Jonas Juselius
    Licensed under terms of LGPLv3 license (see LICENSE-LGPL)
    Repository and issue tracking: https://github.com/rbast/cmake-boot
"""

import subprocess
import re
import os
import sys
from docopt import docopt

__version__ = '0.0.0' # heavy refactoring

#-------------------------------------------------------------------------------

def check_cmake_exists(cmake_command):
    p = subprocess.Popen('%s --version' % cmake_command,
                         shell=True,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE)
    if not ('cmake version' in p.communicate()[0]):
        sys.stderr.write('   This code is built using CMake\n\n')
        sys.stderr.write('   CMake is not found\n')
        sys.stderr.write('   get CMake at http://www.cmake.org/\n')
        sys.stderr.write('   on many clusters CMake is installed\n')
        sys.stderr.write('   but you have to load it first:\n')
        sys.stderr.write('   $ module load cmake\n')
        sys.exit(1)

#-------------------------------------------------------------------------------

def gen_cmake_command(options, arguments):

    command = []

    # take care of compilers
    for lang in ['fc', 'cc', 'cxx']:
        if '--%s' % lang in arguments:
            command.append('%s=%s' % (lang.upper(), arguments['--%s' % lang]))

    command.append('cmake')

    # here we figure out that for instance --omp sets -DENABLE_OMP=ON
    for line in re.findall('.*sets.*ON.*', options):
        flag = line.split()[0]
        action = re.findall('-\w+=ON', line)[0]
        if arguments[flag]:
            command.append(action)

    for libtype in ['blas', 'lapack']:
        arg = arguments['--%s' % libtype]
        if arg == 'builtin':
            command.append('-DENABLE_BUILTIN_%s=ON' % libtype.upper())
            command.append('-DENABLE_AUTO_%s=OFF' % libtype.upper())
        elif arg == 'auto':
            command.append('-DENABLE_AUTO_%s=ON' % libtype.upper())
        elif arg == 'none':
            command.append('-DENABLE_AUTO_%s=OFF' % libtype.upper())
        else:
            if not os.path.isfile(arg):
                sys.stderr.write('ERROR: --%s=%s does not exist\n' % (libtype, arg))
                sys.exit(1)
            command.append('-DEXPLICIT_%s_LIB=%s' % (libtype.upper(), arg))
            command.append('-DENABLE_AUTO_%s=OFF' % libtype.upper())

    # FIXME --mkl and --blas/--lapack conflict
    if arguments['--mkl']:
        possible_mkl_values = ['sequential', 'parallel', 'cluster']
        if arguments['--mkl'] not in possible_mkl_values:
            sys.stderr.write('ERROR: possible --mkl values are: %s\n' % ', '.join(possible_mkl_values))
            sys.exit(1)
        command.append('-DMKL_FLAG="-mkl=%s"' % arguments['--mkl'])
        command.append('-DENABLE_AUTO_BLAS=OFF')
        command.append('-DENABLE_AUTO_LAPACK=OFF')

    if arguments['--explicit-libs']:
        # remove leading and trailing whitespace
        # otherwise CMake complains
        command.append('-DEXPLICIT_LIBS="%s"' % arguments['--explicit-libs'].strip())

    command.append('-DCMAKE_BUILD_TYPE=%s' % arguments['--type'])

    return ' '.join(command)

#-------------------------------------------------------------------------------

def setup_build_path(build_path):
    if os.path.isdir(build_path):
        fname = os.path.join(build_path, 'CMakeCache.txt')
        if os.path.exists(fname):
            sys.stderr.write('aborting setup - build directory %s which contains CMakeCache.txt exists already\n' % build_path)
            sys.stderr.write('remove the build directory and then rerun setup\n')
            sys.exit(1)
    else:
        os.makedirs(build_path, 0755)

#-------------------------------------------------------------------------------

def run_cmake(command, build_path, default_build_path):
    topdir = os.getcwd()
    os.chdir(build_path)
    p = subprocess.Popen(command,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE)
    s = p.communicate()[0]
    # print cmake output to screen
    print(s)
    # write cmake output to file
    f = open('cmake_output', 'w')
    f.write(s)
    f.close()
    # change directory and return
    os.chdir(topdir)
    if 'Configuring incomplete' in s:
        # configuration was not successful
        if (build_path == default_build_path):
            # remove build_path iff not set by the user
            # otherwise removal can be dangerous
            shutil.rmtree(default_build_path)
    else:
        # configuration was successful
        save_configure_command(sys.argv, build_path)
        print_build_help(build_path, default_build_path)

#-------------------------------------------------------------------------------

def print_build_help(build_path, default_build_path):
    print('   configure step is done')
    print('   now you need to compile the sources:')
    if (build_path == default_build_path):
        print('   $ cd build')
    else:
        print('   $ cd ' + build_path)
    print('   $ make')

#-------------------------------------------------------------------------------

def save_configure_command(argv, build_path):
    file_name = os.path.join(build_path, 'configure_command')
    f = open(file_name, 'w')
    f.write(' '.join(argv[:]) + '\n')
    f.close()

#-------------------------------------------------------------------------------

def parse_options(options, argv=None):
    return docopt(options, argv)

#-------------------------------------------------------------------------------

# for testing only
options = """
Usage:
  ./setup [options] [<builddir>]
  ./setup (-h | --help)

Options:
  --cc=<CC>               C compiler [default: gcc].
  --cxx=<CXX>             C++ compiler [default: g++].
  --omp                   Enable OpenMP (sets -DENABLE_OMP=ON).
  --mpi                   Enable MPI (sets -DENABLE_MPI=ON).
  --opencl                Enable OpenCL (sets -DENABLE_OPENCL=ON).
  --coverage              Enable code coverage (sets -DENABLE_CODE_COVERAGE=ON).
  --mkl=<MKL>             Pass MKL flag to the compiler and linker (sequential, parallel, or cluster).
  --blas=<BLAS>           Specify BLAS library (auto, builtin, none, or full path) [default: auto].
  --lapack=<LAPACK>       Specify LAPACK library (auto, builtin, none, or full path) [default: auto].
  --explicit-libs=<LIBS>  Explicit linker specification for extra libraries; passed directly to the linker.
  --type=<TYPE>           Set the CMake build type (debug, release) [default: release].
  <builddir>              Build directory.
  --show                  Show CMake command and exit. Do not write any files.
  -h --help               Show this screen.
"""

#-------------------------------------------------------------------------------

def test_args_to_command_bare():

    arguments = parse_options(options, argv=['./setup'])
    command = gen_cmake_command(options, arguments)
    assert command == 'CC=gcc CXX=g++ cmake -DENABLE_AUTO_BLAS=ON -DENABLE_AUTO_LAPACK=ON -DCMAKE_BUILD_TYPE=release'

#-------------------------------------------------------------------------------

def test_args_to_command():

    arguments = parse_options(options, argv=['./setup', '--cc=icc', '--cxx=icpc', '--mkl=parallel', '--omp'])
    command = gen_cmake_command(options, arguments)
    assert command == 'CC=icc CXX=icpc cmake -DENABLE_OMP=ON -DENABLE_AUTO_BLAS=ON -DENABLE_AUTO_LAPACK=ON -DMKL_FLAG="-mkl=parallel" -DENABLE_AUTO_BLAS=OFF -DENABLE_AUTO_LAPACK=OFF -DCMAKE_BUILD_TYPE=release'

#-------------------------------------------------------------------------------

def test_options_regex():

    arguments = parse_options(options, argv=['./setup', '--omp', '--mpi', '--coverage'])
    command = gen_cmake_command(options, arguments)
    assert command == 'CC=gcc CXX=g++ cmake -DENABLE_OMP=ON -DENABLE_MPI=ON -DENABLE_CODE_COVERAGE=ON -DENABLE_AUTO_BLAS=ON -DENABLE_AUTO_LAPACK=ON -DCMAKE_BUILD_TYPE=release'
