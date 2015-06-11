
# Copyright (c) 2015 by Radovan Bast and Jonas Juselius
# See https://github.com/scisoft/autocmake/blob/master/LICENSE


import subprocess
import os
import sys
import shutil


def check_cmake_exists(cmake_command):
    """
    Check whether CMake is installed. If not, print
    informative error message and quits.
    """
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


def setup_build_path(build_path):
    """
    Create build directory. If this already exists, print informative
    error message and quit.
    """
    if os.path.isdir(build_path):
        fname = os.path.join(build_path, 'CMakeCache.txt')
        if os.path.exists(fname):
            sys.stderr.write('aborting setup\n')
            sys.stderr.write('build directory %s which contains CMakeCache.txt already exists\n' % build_path)
            sys.stderr.write('remove the build directory and then rerun setup\n')
            sys.exit(1)
    else:
        os.makedirs(build_path, 0755)


def run_cmake(command, build_path, default_build_path):
    """
    Execute CMake command.
    """
    topdir = os.getcwd()
    os.chdir(build_path)
    p = subprocess.Popen(
            command,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE
        )
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


def print_build_help(build_path, default_build_path):
    """
    Print help text after configuration step is done.
    """
    print('   configure step is done')
    print('   now you need to compile the sources:')
    if (build_path == default_build_path):
        print('   $ cd build')
    else:
        print('   $ cd ' + build_path)
    print('   $ make')


def save_configure_command(argv, build_path):
    """
    Save configure command to a file.
    """
    file_name = os.path.join(build_path, 'configure_command')
    f = open(file_name, 'w')
    f.write(' '.join(argv[:]) + '\n')
    f.close()


def configure(root_directory, build_path, cmake_command, only_show):
    """
    Main configure function.
    """
    default_build_path = os.path.join(root_directory, 'build')

    # check that CMake is available, if not stop
    check_cmake_exists('cmake')

    # deal with build path
    if build_path is None:
        build_path = default_build_path
    if not only_show:
        setup_build_path(build_path)

    print('%s\n' % cmake_command)
    if only_show:
        sys.exit(0)

    run_cmake(cmake_command, build_path, default_build_path)