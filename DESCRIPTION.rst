binstarsolver
=============

Overview
^^^^^^^^

The purpose of this project is to estimate physical quantities of an eclipsing, spectroscopic binary star system using observed quantities from photometry and spectroscopy.
The estimated physical quantities serve as input to a more complete numerical model such as that of Wilson-Divinney, c.f. chapter 7 of [3]_.
The model used in this project is a simple, spherical detached binary star model with no orbital eccentricity and no limb darkening.
The model is from chapter 7 of [1]_ and chapter 7 of [2]_.

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

References
----------

.. [1] `Budding, 2007, Introduction to Astronomical Photometry <https://books.google.com/books?id=g_K3-bQ8lTUC>`_
.. [2] `Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics <https://books.google.com/books?id=M8wPAQAAMAAJ>`_
.. [3] `Kallrath and Milone, 2007, Eclipsing Binary Stars: Modeling and Analysis <https://books.google.com/books?id=CrXBnZFdjXgC>`_
