[![Build Status](https://travis-ci.org/dftlibs/xcint.svg?branch=master)](https://travis-ci.org/dftlibs/xcint/builds)
[![Coverage Status](https://coveralls.io/repos/dftlibs/xcint/badge.png?branch=master)](https://coveralls.io/r/dftlibs/xcint?branch=master)
[![Documentation Status](https://readthedocs.org/projects/xcint/badge/?version=latest)](http://xcint.readthedocs.org)


XCint
=====

Exchange-correlation integrator.

- [Documentation](http://xcint.readthedocs.org/)
- [Build and test history](https://travis-ci.org/dftlibs/xcint/builds)
- [Code coverage](https://coveralls.io/r/dftlibs/xcint)
- Built with [Autocmake](https://github.com/scisoft/autocmake)
- Status: pre-alpha
- Licensed under [BSD-3](../master/LICENSE)

Primary test environments
=========================

Continuous integration builds
-----------------------------

The Travis CI builds are triggered only when pushing to the `master` branch.
All Travis CI builds on master use ccache to speed up execution.

- Ubuntu 12.04 LTS 64-bit with Python 2.7.3 and CMake 3.3.2
  this is the environment offered by [Travis CI](https://travis-ci.org) pulling
  in various PPA. The following compilers are used, both in release and debug:

  1. GCC 4.6
  2. GCC 4.7
  3. GCC 4.8
  4. GCC 4.9
  5. GCC 5.1, with and without coverage analysis in debug mode.
     Coverage analysis performed using [Coveralls](https://coveralls.io/)
  6. Clang 3.5 and GFortran 4.6
  7. Clang 3.6 and GFortran 4.6
  8. Clang 3.7 and GFortran 4.6
  9. Clang 3.8 and GFortran 4.6

- Mac OS X 10.9.5 with Python 2.7.10 and CMake 3.2.3
  this is the environment offered by [Travis CI](https://travis-ci.org)
  The following compilers are used, both in release and debug:

  1. XCode 6.4 with Clang and GFortran 5.2
  2. XCode 6.4 with GCC 5.2
  3. XCode 7.0 with Clang and GFortran 5.2
  4. XCode 7.0 with GCC 5.2
