binstarsolver
=============

.. image:: https://zenodo.org/badge/4765/ccd-utexas/binstarsolver.svg   
        :target: http://dx.doi.org/10.5281/zenodo.16983

Overview
--------

The purpose of this project is to estimate physical quantities of an eclipsing, spectroscopic binary star system using observed quantities from photometry and spectroscopy.
The estimated physical quantities serve as input to a more complete numerical model such as that of Wilson-Devinney, c.f. chapter 7 of [1]_.
The model used in this project is a simple, spherical detached binary star model with no orbital eccentricity and no limb darkening.
The model is from chapter 7 of [2]_ and chapter 7 of [3]_.

Input observed quantities
-------------------------

* Relative times of begin/end ingress/egress and mid-eclipse.
* Orbital period.
* Light levels during occultation (primary eclipse), transit (secondary eclipse), and outside of eclipse.
* Stellar radial velocities.

Output estimated physical quantities
------------------------------------

* Orbital inclination angle.
* Orbital radii.
* Stellar masses.
* Stellar radii.
* Ratio of stellar radiative fluxes.
* Ratio of stellar effective temperatures.
* Ratio of stellar luminosites.

Examples and Q&A forum
----------------------

* For explanations of examples, please see the `binstarsolver GitHub wiki <https://github.com/ccd-utexas/binstarsolver/wiki>`_.
* For questions and discussion, please see the `binstarsolver Google group <https://groups.google.com/forum/#!forum/binstarsolver>`_.
* Got a second for a `user survey <https://docs.google.com/forms/d/1vneANTMMaOdQSRvIm2OJYItgaTTbp4f9EM8ImKqwD-g/viewform>`_? I'd love your feedback! Thank you for being part of open science. (All survey responses are anonymous.)

Citations
---------

* If this code is useful to your academic research, please consider citing it by including this DOI link: http://dx.doi.org/10.5281/zenodo.16983

Recent changes
--------------

* [v0.1.3] - 20150420T180000
   - Created Google group, DOI for citations, CHANGELOG.rst, IPython Notebook examples.
   - Fixed special cases where solver for inclination did not converge.
* [v0.1.2] - 20150123T172000
   - Initial release to PyPI.

References
----------

.. [1] `Kallrath and Milone, 2007, Eclipsing Binary Stars: Modeling and Analysis <https://books.google.com/books?id=CrXBnZFdjXgC>`_
.. [2] `Budding, 2007, Introduction to Astronomical Photometry <https://books.google.com/books?id=g_K3-bQ8lTUC>`_
.. [3] `Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics <https://books.google.com/books?id=M8wPAQAAMAAJ>`_
