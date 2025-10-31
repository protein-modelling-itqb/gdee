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

Core Installation
-----------------

1. **Create Python Virtual Environment**

.. code-block:: bash

    # Create and activate virtual environment
    python3 -m venv gdee-env
    source gdee-env/bin/activate

2. **Install GDEE Platform**

.. code-block:: bash

    # install from source
    git clone https://github.com/protein-modelling-itqb/gdee.git
    cd gdee/package
    pip install -e .

3. **Verify Core Installation**

.. code-block:: python

    # Test basic import
    import gdee
    from gdee import ProteinEngineering
    print("GDEE platform installed successfully")

Required External Programs
--------------------------

The GDEE platform requires several external programs for structure modeling and molecular docking.

MODELLER Installation
~~~~~~~~~~~~~~~~~~~~~

MODELLER is required for 3D structure modeling of protein variants.

1. **Register and Download**
   
   - Visit https://salilab.org/modeller/
   - Register for academic license
   - Download MODELLER for your platform

2. **Install MODELLER**

.. code-block:: bash

    # Ubuntu/Debian
    sudo apt-get install modeller

    # Or from downloaded package
    sudo dpkg -i modeller_*.deb

3. **Configure License**

.. code-block:: bash

    # Edit MODELLER configuration
    sudo nano /usr/lib/modeller*/modlib/modeller/config.py
    
    # Add your license key
    license = 'YOUR_LICENSE_KEY'

4. **Install Python Bindings**

.. code-block:: bash

    # In your GDEE virtual environment
    pip install modeller

MGLTools Installation
~~~~~~~~~~~~~~~~~~~~~

MGLTools is required for converting PDB files to PDBQT format for docking.

1. **Download MGLTools**

.. code-block:: bash

    # Download MGLTools 1.5.6
    wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz
    
    # Extract
    tar -xzf mgltools_x86_64Linux2_1.5.6.tar.gz
    
    # Install
    cd mgltools_x86_64Linux2_1.5.6
    ./install.sh

2. **Set Installation Path**

.. code-block:: python

    # In your GDEE script
    eng.programs["mgltools"] = "/path/to/mgltools_x86_64Linux2_1.5.6"

AutoDock Vina Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~

AutoDock Vina is required for molecular docking calculations.

1. **Install from Package Manager**

.. code-block:: bash

    # Ubuntu/Debian
    sudo apt-get install autodock-vina


2. **Or Compile from Source**

.. code-block:: bash

    # Download and compile
    wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.0/vina_1.2.0_linux_x86_64
    chmod +x vina_1.2.0_linux_x86_64
    sudo mv vina_1.2.0_linux_x86_64 /usr/local/bin/vina

3. **Configure Path**

.. code-block:: python

    # In your GDEE script
    eng.programs["vina"] = "/usr/local/bin/vina"

Vinardo Installation (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vinardo provides improved docking accuracy compared to standard Vina.

1. **Download Vinardo**

.. code-block:: bash

    # Download from official source
    wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.0/vinardo_1.2.0_linux_x86_64
    chmod +x vinardo_1.2.0_linux_x86_64
    sudo mv vinardo_1.2.0_linux_x86_64 /usr/local/bin/vinardo

2. **Configure Path**

.. code-block:: python

    # In your GDEE script
    eng.programs["vinardo"] = "/usr/local/bin/vinardo"

VoroMQA Installation
~~~~~~~~~~~~~~~~~~~~

VoroMQA is used for model quality assessment using Voronoi tessellation.

1. **Install Voronota Package**

.. code-block:: bash

    # Download Voronota
    wget https://github.com/kliment-olechnovic/voronota/releases/download/v1.25.3049/voronota_1.25.3049.tar.gz
    tar -xzf voronota_1.25.3049.tar.gz
    cd voronota_1.25.3049
    
    # Compile
    make
    
    # Install
    sudo cp voronota-voromqa /usr/local/bin/

2. **Configure Path**

.. code-block:: python

    # In your GDEE script
    eng.programs["voromqa"] = "/usr/local/bin/voronota-voromqa"


Next Steps
----------

After successful installation:

1. Review the :doc:`usage_instructions` for workflow examples
2. Check the :doc:`configuration_reference` for detailed parameter descriptions
3. See the :doc:`troubleshooting` guide for common issues
4. Explore the :doc:`api_reference` for advanced usage