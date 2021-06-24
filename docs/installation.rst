Installation
============

HAYSTAC can be run on either macOS or Linux based systems, and the source code is available on ``github``.

The recommended way to install ``haystac``, and all if its dependencies, is via the [mamba](https://mamba.readthedocs.io/en/latest/installation.html)
package manager (a fast replacement for [conda](https://docs.conda.io/projects/conda/en/latest/index.html)).

You can install ``haystac`` using ``conda``, however, it will take significantly longer to install and analyses will run slower.

Install mamba
---------------------
If you do not have either ``mamba`` or ``conda`` already installed, please refer to the [install instructions](
https://mamba.readthedocs.io/en/latest/installation.html) for ``mambaforge``.

If you have ``conda`` installed, but not ``mamba``, then install ``mamba`` into the base environment::

    conda install -n base -c conda-forge mamba

Install haystac
---------------------
Then use ``mamba`` to install ``haystac`` into a new environment::

    mamba create -c bioconda -n haystac haystac

And activate the environment::

    conda activate haystac

We recommend that you install ``haystac`` into a new environment to avoid dependency conflicts with other software.

Git
---------------------

Clone from ``github``::

    git clone https://github.com/antonisdim/haystac.git


