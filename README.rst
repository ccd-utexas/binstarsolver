binstarsolver
=============

Estimate physical quantities of a binary star system from observed quantities.

Description
-----------

See DESCRIPTION.rst.

.. toctree::

   DESCRIPTION

Examples 
--------

For explanations of examples, please see the `wiki <https://github.com/ccd-utexas/binstarsolver/wiki>`_.
Contributions are always welcome!

Questions
---------

For questions and discussion, please see the `Google group <https://groups.google.com/forum/#!forum/binstarsolver>`_.
Contributions are always welcome!

Installation
------------

Use `Python Packaging Index <https://pypi.python.org/pypi>`_:

::

   $ pip install binstarsolver

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

See CONTRIBUTING.rst.

.. toctree::

   CONTRIBUTING
