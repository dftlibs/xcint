

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


Fortran interface
-----------------

An optional Fortran interface is currently broken and will be restored/reintroduced
as soon as the C interface stabilizes a bit.


Outlook
-------

The plan is to also move the density evaluation (together with AO evaluation
and "Fock"-type matrix distribution routines) outside of XCint to a separate
library.  This will make XCint very thin and compact.  Another advantage is
that this step will make XCint basis-set agnostic and point-group symmetry
agnostic and therefore more general.
