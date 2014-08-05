

[![Build Status](https://travis-ci.org/rbast/cmake-boot.svg?branch=master)](https://travis-ci.org/rbast/cmake-boot/builds)

CMake-boot
==========

Bootstrapping front-end to CMake.

- Uses [docopt](http://docopt.org).
- Copyright 2014 by Radovan Bast and Jonas Jusélius.
- Status: Pre-alpha. Both interfaces and internals subject to heavy modifications.
- Licensed under [LGPLv3](../master/LICENSE).


## Example

First clone the code:
```
$ git clone https://github.com/rbast/cmake-boot.git cboot
```

Then:
```python
#!/usr/bin/env python

import os
from cboot.cboot import configure

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
  --show                  Show CMake command and exit. Do not write any files.
  <builddir>              Build directory.
  -h --help               Show this screen.
"""

root_directory = os.path.dirname(os.path.realpath(__file__))

configure(options, root_directory)
```
