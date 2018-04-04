

What XCint does
===============

XCint integrates the exchange-correlation (XC) energy :math:`E_\text{xc}`
and the elements of the XC potential matrix :math:`V_\text{xc}`,
as well as their derivatives with respect to electric field and/or geometric
perturbations. The integration is performed on a standard numerical grid.

Before the XC energy can be computed, we require the densities.
This is done in two steps:

.. math::

  n_b = \sum_k \chi_{kb} \sum_l D_{kl} \chi_{lb} = \sum_k \chi_{kb} X_{kb}

.. math::

  X_{kb} = \sum_l D_{kl} \chi_{lb}

Then the XC energy is computed using the XC energy density :math:`\epsilon_\text{xc}`
evaluated with the help of XCFun:

.. math::

  E_\text{xc} = \sum_b w_b \epsilon_\text{xc} (n_b)

A similar strategy is used to compute :math:`V_\text{xc}` matrix elements, again in two steps:

.. math::

  (V_\text{xc})_{kl} = \sum_b w_b \chi_{kb} v_\text{xc} (n_b) \chi_{lb}
                     = \sum_b     W_{kb}                      \chi_{lb}

.. math::

  W_{kb} = w_b \chi_{kb} v_\text{xc} (n_b)

Note that XCFun is called only once and returns :math:`\epsilon_\text{xc}`
and :math:`v_\text{xc}` in one go.

In the above scheme we work with batches of points in order to
exploit vectorization and screening, making use of BLAS level 3 libraries
to compute :math:`X_{kb}` and :math:`(V_\text{xc})_{kl}`.
