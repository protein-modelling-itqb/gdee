========================================
  Gene Discovery and Enzyme Engineering
========================================

.. figure:: docs/source/_static/images/cover2.png
   :alt: GDEE Platform Banner
   :width: 800px
   :align: center

**GDEE** is a comprehensive Python package for protein engineering workflows that provides the functionality necessary to run the Gene Discovery and Enzyme Engineering Platform developed under the project ShikiFactory 100 (European Union's Horizon 2020 research and innovation programme under grant agreement number 814408).

GDEE automates the complete computational pipeline for testing several protein variants, from sequence generation to molecular docking evaluation and rescoring, making it ideal for enzyme engineering and gene discovery projects.

Features
========

**Variant Generation**
- **MSA-based variants**: Generate variants from FASTA sequences (BLAST results)
- **Random mutations**: BLOSUM62-guided amino acid substitutions
- **Exhaustive mutations**: Systematic exploration of mutation combinations
- **Flexible selection**: Target specific residues with exclusion rules

 **3D Structure Modeling**
- **MODELLER integration**: Homology modeling with mutation-specific optimization
- **Quality assessment**: VoroMQA and Normalized DOPE scoring
- **Batch processing**: Multiple models per variant with automatic selection

 **Molecular Docking**
- **AutoDock Vina/Vinardo**: High-throughput protein-ligand docking
- **Multiple ligands**: Parallel evaluation of different compounds
- **Pose analysis**: Energy ranking and geometric measurements

 **Analysis & Storage**
- **SQLite database**: Comprehensive result storage with relationships
- **Filtering**: SQL-based result analysis and ranking
- **Export capabilities**: CSV and database export for further analysis
- **Measurement System**: Distance calculations and geometric analysis

 **Scalable Execution**
- **Single machine**: Simple sequential processing
- **MPI support**: Distributed computing across clusters
- **File archiving**: Automatic result compression and organization

 **Rescoring**
- **Post-hoc Analysis**: Re-score existing docking results with trained metamodel
- **Result Ranking**: Quickly identify top candidates using alternative scoring

How to Install GDEE
===================

GDEE requires several external dependencies including MODELLER, AutoDock Vina, and MGLTools for complete functionality. The installation process involves setting up both the Python package and external computational tools.

**System Requirements**
- Python 3.6+
- Linux
- MODELLER license (academic license available)

For detailed installation instructions including external program setup, and platform-specific instructions, please see the complete **Installation Guide** in our documentation:

 **Installation Documentation**: https://gdee.readthedocs.io/en/latest/installation.html

How to Use GDEE
===============

GDEE provides a high-level Python interface for protein engineering workflows. Here's a basic example to get you started:

**Basic Workflow Example**

.. code-block:: python

    from gdee import ProteinEngineering

    # Initialize workflow
    eng = ProteinEngineering("MyProtein", "results.db")

    # Configure template structure
    eng.pdb = "template.pdb"

    # Set up variant generation
    eng.variant["name"] = "mutation"
    eng.variant["selection"] = "A:123 A:456 B:789"  # Target specific residues
    eng.variant["max_iterations"] = 100
    eng.variant["combinations"] = 2

    # Configure 3D modeling
    eng.model["num_models"] = 5
    eng.model["optimize_level"] = 1

    # Add ligand for docking
    ligand = eng.add_ligand("compound1", "ligand.pdbqt")

    # Set docking parameters
    eng.evaluator["name"] = "vina"
    eng.evaluator["box_center"] = [10.0, 15.0, 20.0]
    eng.evaluator["box_size"] = [20.0, 20.0, 20.0]

    # Run the complete workflow
    eng.run()

**Rescoring Workflow Example**

Re-evaluate existing docking results with trained metamodel:

.. code-block:: python

    from gdee import RescoreVariants, Ranked_Database

    # Open existing results database as a Ranked_Database-like object
    ranked_db = Ranked_Database("results.db")

    # List of trained scoring function model files
    scoring_functions = ["scoring_functions.pkl"]

    # Perform rescoring on existing results
    rescorer = RescoreVariants(
        in_db=ranked_db,
        out_db="rescored_results.db",
        table="rescored_poses",
        functions=scoring_functions,
        files_path="/path/to/structure/files"
    )

    rescorer.run()

For comprehensive tutorials, advanced configuration options, workflow examples, and best practices, please refer to our complete **Usage Documentation**:

 **Usage Guide**: https://gdee.readthedocs.io/en/latest/usage.html


How Does GDEE Work
==================

GDEE implements a modular pipeline architecture that processes protein variants through multiple computational stages in an automated workflow. The platform integrates variant generation, structure modeling, quality assessment, molecular docking, and analysis into a unified framework.

**Pipeline Architecture**

The GDEE workflow consists of the following stages:

1. **Variant Generation**: Creates protein variants using a FASTA file, random mutations, or exhaustive mutations
2. **Structure Modeling**: Uses MODELLER for homology modeling with mutation-specific optimization
3. **Quality Assessment**: Evaluates model quality using VoroMQA and Normalized DOPE scores
4. **Molecular Docking**: Performs protein-ligand docking with AutoDock Vina/Vinardo
5. **Measurement Calculation**: Computes distances and geometric metrics
6. **Results Storage**: Saves all data to SQLite database for analysis

**Rescoring Module**

The rescoring feature enables efficient re-evaluation of existing docking poses:

- Reads poses from the ranked database
- Applies trained metamodel without structural changes
- Stores new scores in a database for comparison

For detailed information about the algorithmic approaches, computational methods, validation studies, and scientific applications of the GDEE platform, please refer to our research publication:

 **Scientific Publication**: https://doi.org/10.1101/2025.09.09.675117

Quick Links
===========

 **Complete Documentation**: https://gdee.readthedocs.io/

- `Installation Guide <https://gdee.readthedocs.io/en/latest/installation.html>`_
- `Usage Tutorials <https://gdee.readthedocs.io/en/latest/usage.html>`_
- `Configuration Reference <https://gdee.readthedocs.io/en/latest/configuration.html>`_

 **Scientific Publication**: https://doi.org/10.1101/2025.09.09.675117


Citing GDEE
===========

If you use GDEE in your research, please cite:

.. code-block:: bibtex

    @article{souza2025gdee,
      title={GDEE: A Structure-Based Platform for Gene Discovery and Enzyme Engineering},
      author={Souza, Caio S and Correia, Jo{\~a}o PG and Rocha, Isabel and Lousa, Diana and Soares, Cl{\'a}udio M},
      journal={bioRxiv},
      pages={2025--09},
      year={2025},
      publisher={}
    }

Contributing
============

We welcome contributions! Please see our documentation for guidelines on:

- Bug reports and feature requests
- Development setup and testing
- Code contribution workflow
- Documentation improvements

License
=======

.. |cc-by-nc-sa-shield| image:: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
   :target: http://creativecommons.org/licenses/by-nc-sa/4.0/
   :alt: CC BY-NC-SA 4.0

.. |cc-by-nc-sa-image| image:: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
   :target: http://creativecommons.org/licenses/by-nc-sa/4.0/
   :alt: CC BY-NC-SA 4.0

|cc-by-nc-sa-shield|

This work is licensed under a
`Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License`_.

|cc-by-nc-sa-image|

.. _Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License: http://creativecommons.org/licenses/by-nc-sa/4.0/

For commercial use or licensing inquiries, please contact:

- **Corresponding Author**: Cláudio M. Soares (claudio@itqb.unl.pt)
- **Institution**: Instituto de Tecnologia Química e Biológica António Xavier

Acknowledgments
===============

This work was supported by the ShikiFactory 100 project under the European Union's Horizon 2020 research and innovation programme (grant agreement number 814408).
