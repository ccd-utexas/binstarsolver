binstarsolver
=============

Estimate physical quantities of an eclipsing, spectroscopic binary star system from observed lightcurves and radial velocities.

Description
-----------

Overview
^^^^^^^^

This model is a simple spherical detached binary star model with no orbital eccentricity, no limb darkening. The model is from chapter 7 of [1]_ and chapter 7 of [2]_. The purpose of this package is to generate initial parameters as input to a more complete numerical model such as that of Wilson-Divinney, c.f. chapter 7 of [3]_.

Input observable quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Relative times of begin/end ingress/egress and mid-eclipse.
* Orbital period.
* Light levels during occultation (primary eclipse), transit (secondary eclipse), outside of eclipse.
* Stellar radial velocities.

Output estimated physical quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Orbital inclination angle.
* Orbital radii.
* Stellar masses.
* Stellar radii.
* Ratio of stellar radiative fluxes.
* Ratio of stellar effective temperatures.
* Ratio of stellar luminosites.

Examples
--------

* See repository `wiki <https://github.com/ccd-utexas/binstarsolver/wiki>`_.

Installation
------------

See INSTALL.md.

Contributing
------------

See CONTRIBUTING.md.

References
----------

.. [1] `Budding, 2007, Introduction to Astronomical Photometry <https://books.google.com/books?id=g_K3-bQ8lTUC>`_
.. [2] `Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics <https://books.google.com/books?id=M8wPAQAAMAAJ>`_
.. [3] `Kallrath and Milone, 2007, Eclipsing Binary Stars: Modeling and Analysis <https://books.google.com/books?id=CrXBnZFdjXgC>`_
