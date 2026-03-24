Installation and Setup
======================

This guide provides step-by-step instructions for installing the GDEE platform and all required dependencies.

System Requirements
-------------------

**Operating System:**
- Linux (Ubuntu 18.04+, or equivalent)

**Python Environment:**
- Python 3.7 or higher
- pip package manager
- Virtual environment support (recommended)

**Python Package Requirements**
- numpy>=1.14
- mdanalysis>=0.20
- mpi4py>=3.0
- path>=13.0
- biopython>=1.75
- oddt>=0.7
- openbabel-wheel>=3.1.1.1
- six >= 1.17.0
- scipy == 1.9.0

Core Installation
-----------------

1. **Create Python Virtual Environment**

.. code-block:: bash

    # Create and activate virtual environment
    python3 -m venv gdee-env
    source gdee-env/bin/activate

2. **Install six and scipy**

    .. note::
        scipy==1.9.0 is required for compatibility with the oddt package

.. code-block:: bash

    pip install six>=1.17.0 scipy==1.9.0

3. **Install GDEE Platform**

.. code-block:: bash

    # install from source
    git clone https://github.com/protein-modelling-itqb/gdee.git
    cd gdee/package
    pip install -e .


Required External Programs
--------------------------

The GDEE platform requires several external programs. Please follow the official installation documentation for each tool to ensure you have the most up-to-date installation procedures.

MODELLER Installation
~~~~~~~~~~~~~~~~~~~~~

MODELLER is required for 3D structure modeling of protein variants.

**Installation:** Please follow the official MODELLER installation guide at https://salilab.org/modeller/download_installation.html

**Requirements:**
- Academic license registration required
- Python bindings must be installed in your GDEE virtual environment


MGLTools Installation
~~~~~~~~~~~~~~~~~~~~~

MGLTools is required for converting PDB files to PDBQT format for docking.

**Installation:** Please follow the official MGLTools installation guide at https://ccsb.scripps.edu/mgltools/downloads/ 

**Configuration in GDEE:**

.. code-block:: python

    # In your GDEE script, specify the MGLTools installation path
    eng.programs["mgltools"] = "/path/to/mgltools_installation"

AutoDock Vina Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~

AutoDock Vina is required for molecular docking calculations.

**Installation:** Please follow the official AutoDock Vina installation guide at https://github.com/ccsb-scripps/AutoDock-Vina

**Configuration in GDEE:**

.. code-block:: python

    # In your GDEE script, specify the Vina executable path
    eng.programs["vina"] = "/path/to/vina"

Vinardo Installation (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vinardo is available as a scoring function in the Smina program.

**Installation:** Please follow the official Smina installation guide at https://sourceforge.net/projects/smina/

**Configuration in GDEE:**

.. code-block:: python

    # In your GDEE script, specify the Smina executable path
    eng.programs["vinardo"] = "/path/to/smina"

VoroMQA Installation (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VoroMQA is used for model quality assessment using Voronoi tessellation.

**Installation:** Please follow the official Voronota installation guide at https://github.com/kliment-olechnovic/voronota

**Configuration in GDEE:**

.. code-block:: python

    # In your GDEE script, specify the VoroMQA executable path
    eng.programs["voromqa"] = "/path/to/voronota-voromqa"


Next Steps
----------

After successful installation:

1. Review the :doc:`Usage` for workflow examples
2. Ensure all external programs are properly configured in your GDEE scripts
3. Test your installation with a simple workflow to verify all dependencies work correctly