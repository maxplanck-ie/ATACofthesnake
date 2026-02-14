Installation
============

Setup
-----

While not strictly necessary, it is recommended to set up a `Snakemake profile <https://snakemake.readthedocs.io/en/stable/executing/cli.html#executing-profiles>`_
to run the workflow. This way submission via a scheduler is possible, which (drastically) reduces the runtime versus a local run.
Additionally, the environments within the workflow are managed via `Conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.
When using a Snakemake profile, make sure the `use-conda` option is set to `True`. If you are not using a profile, Snakemake is called with the `use-conda` option automatically.
Conda needs to be installed and available in your PATH, regardless of the way you run the workflow (either profile-based or not).

Recommended
-----------

ATACofthesnake is available on PyPI, and installation is recommended via pip:

.. code:: bash

   pip install ATACofthesnake

or via uv

.. code:: bash

   uv pip install ATACofthesnake

Development version
-------------------

You can install the latest development version directly via GitHub:

.. code:: bash

    pip install git+https://github.com/maxplanck-ie/ATACofthesnake.git

or again via uv:

.. code:: bash

    uv pip install git+https://github.com/maxplanck-ie/ATACofthesnake.git
