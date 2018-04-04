

Status and outlook
==================


Status
------

The API is pre-alpha. Expect significant changes.


Grid generation
---------------

The grid generation has been moved outside XCint
(the tests employ https://github.com/dftlibs/numgrid).
This has the following advantages:

- Gives the caller the possibility to use other grid generators.
- Makes grid-based (MPI) parallelization relatively trivial.
- Moves grid-based (MPI) parallelization outside XCint.


MPI parallelization
-------------------

MPI parallelization has been removed (see above section) as it can be
introduced by the caller in very few lines.  Not having MPI parallelization
inside XCint simplifies the code and testing.


Functional parsing
------------------

This is currently done inside the code but should move outside. The reason why the functional
is parsed and tracked inside the code is that XCFun does not allow to track the same functional
using both the Fortan and C interfaces in the same run. Moving the functional parsing out now
would break the Fortran interface of XCint.


Outlook
-------

The plan is to also move the density evaluation (together with AO evaluation
and "Fock"-type matrix distribution routines) outside of XCint to a separate
library.  This will make XCint very thin and compact.  Another advantage is
that this step will make XCint basis-set agnostic and point-group symmetry
agnostic and therefore more general.
