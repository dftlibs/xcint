

Interfacing
===========


General hints
-------------

The first is that the density matrix has to be scaled by half before calling
the routine, and that is already counting that D is spin-restricted.
Similarly, the output F needs to be multiplied times two. This is probably
historically related to how DALTON handles data.

The AO primitive weights are expected to contain the spherical harmonic
normalization constant. This is essential because there are different ways of
absorbing the angular and radial normalizations into different quantities,
depending on definitions.
