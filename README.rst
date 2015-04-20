binstarsolver
=============

.. image:: https://zenodo.org/badge/4765/ccd-utexas/binstarsolver.svg   
        :target: http://dx.doi.org/10.5281/zenodo.16983

Estimate physical quantities of a binary star system from observed quantities.

Description
-----------

For the Python Packaging Index `binstarsolver <https://pypi.python.org/pypi/binstarsolver>`_ description, see DESCRIPTION.rst.

.. toctree::

   DESCRIPTION

Examples and Q&A forum
----------------------

* For explanations of examples, please see the `binstarsolver GitHub wiki <https://github.com/ccd-utexas/binstarsolver/wiki>`_.
* For questions and discussion, please see the `binstarsolver Google group <https://groups.google.com/forum/#!forum/binstarsolver>`_.
* Got a second for a `user survey <https://docs.google.com/forms/d/1vneANTMMaOdQSRvIm2OJYItgaTTbp4f9EM8ImKqwD-g/viewform>`_? I'd love your\
 feedback! Thank you for being part of open science. (All survey responses are anonymous.)

Citations
---------

* If this code is useful to your academic research, please consider citing it by including this DOI link: http://dx.doi.org/10.5281/zenodo.16983

Installation and upgrades
-------------------------

Use `Python Packaging Index <https://pypi.python.org/pypi>`_:

::

   $ pip install binstarsolver
   $ pip install binstarsolver --upgrade

Testing
-------

Use `pytest <https://pytest.org>`_:

::

   $ git clone https://github.com/ccd-utexas/binstarsolver.git
   $ cd binstarsolver
   $ git tag --list
   $ git checkout tags/v0.1.2 # or latest tag name
   $ pip install -e .[test] # installs packages for testing
   $ py.test -v

Contributing
------------

For contributing code, please see CONTRIBUTING.rst.
For other ways to contribute, please see the sections 'Examples' and 'Q&A forum' above.

.. toctree::

   CONTRIBUTING
