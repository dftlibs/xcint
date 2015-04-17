

Interface to XCFun
==================


Functional definition
---------------------

The functional class is defined in ``src/Functional.h``.
The functional is set using the ``set_functional`` function.

The following example sets the LDA functional:

.. code-block:: c

  double hfx, mu, beta;
  set_functional("lda", hfx, mu, beta);

Here ``hfx``, ``mu``, and ``beta`` are set by the ``set_functional`` function.

The order of the differentiation must be set using the function ``set_order``:

.. code-block:: c

  set_order(1); // in this case first-order differentiation


Evaluating functional values and derivatives
--------------------------------------------

The function ``xc_eval`` receives a functional object (``fun.fun``)
and a Taylor series expansion of the density variables (``xcin``)
and computes a Taylor series expansion of the XC energy density.
