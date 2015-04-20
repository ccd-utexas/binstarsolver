CONTRIBUTING
============

Thanks for your interest in contributing! Please use the following:

* GitHub guidelines: https://guides.github.com/
* Git branching model: http://nvie.com/posts/a-successful-git-branching-model/
* Docstring style: https://github.com/numpy/numpy/blob/master/doc/example.py
* Coding style: https://google-styleguide.googlecode.com/svn/trunk/pyguide.html
* Versioning: http://semver.org/
* Change log: http://keepachangelog.com/
* Packaging: https://packaging.python.org/en/latest/distributing.html

How to issue a release
----------------------
* Edit the project in develop mode: ``$ python setup.py develop``
* Complete the issues for the release milestone, and close the milestone: https://github.com/ccd-utexas/binstarsolver/milestones
* Edits should be merged with ``develop`` branch.
* Run tests with ``$ py.test -v`` and run the supported examples from ``binstarsolver/examples``
* TODO: Check the travis ci build https://travis-ci.org/ccd-utexas/binstarsolver
* Create a branch ``release_vx.x.x`` off of ``develop`` branch.
* Run tests again with ``$ py.test -v`` and run the supported examples from ``binstarsolver/examples``
* Edits should be merged with ``release_vx.x.x`` branch.
* Merge ``release_vx.x.x`` branch into ``master``.
* TODO: check the travis ci build again https://travis-ci.org/ccd-utexas/binstarsolver
* Edits should be merged with ``master`` branch.
* Create a test distribution: ``$ python setup.py sdist upload -r https://testpypi.python.org/pypi``
* TODO: make a wheel distribution
* TODO: use twine
* Inspect website and links: https://testpypi.python.org/pypi/binstarsolver
* Test installation: ``$ pip uninstall binstarsolver ; pip install -i https://testpypi.python.org/pypi binstarsolver``
* Create a release vx.x.x on GitHub and name with a date-timestamp: https://github.com/ccd-utexas/binstarsolver/releases
* Update the DOI badge and links from zenodo.org in DESCRIPTION.rst, README.rst
* Create release distribution: ``$ python setup.py sdist upload``
* Merge ``master`` branch into ``release_vx.x.x``
* Merge ``release_vx.x.x`` into ``develop``
* Delete ``release_vx.x.x`` branch
