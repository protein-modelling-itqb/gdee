Usage Instructions
==================

This section provides practical examples demonstrating how to run the Gene Discovery and Enzyme Engineering features of the GDEE Platform.

Overview
--------

The GDEE platform provides a unified interface through the :class:`~gdee.engineer.ProteinEngineering` class for executing comprehensive protein engineering workflows. The platform supports:

- Multiple variant generation strategies (from sequence-based gene discovery to mutation-based enzyme engineering)
- 3D structure modeling with MODELLER
- Model quality assessment using VoroMQA and Normalized DOPE
- Molecular docking with AutoDock Vina/Vinardo
- Distance measurements for pose filtering
- Both single-machine and distributed MPI execution
- SQLite database storage for all results

The following sections illustrate how to configure and run typical workflows for gene discovery, enzyme engineering, and rescoring existing results.


Configuration Parameters Reference
----------------------------------

This section provides a comprehensive reference for all configuration parameters available in the GDEE platform through the :class:`~gdee.engineer.ProteinEngineering` class.

Platform Configuration
~~~~~~~~~~~~~~~~~~~~~~

Control the execution environment and computational resources:

.. code-block:: python

    # Platform execution settings
    eng.platform["name"] = "simple"     # Execution mode: "simple" or "mpi"
    eng.platform["local_cpu"] = 1       # Number of CPU cores per MPI process

**Available Options:**

- **name**:
  - ``"simple"``: Single-threaded execution on local machine
  - ``"mpi"``: Distributed execution using MPI for parallel processing
- **local_cpu**: Integer specifying CPU cores for local parallelism within each MPI process

File Management Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure output file organization and compression:

.. code-block:: python

    # File archiving and output management
    eng.io["output"] = "files"                   # Base directory for output files
    eng.io["output_format"] = ".{:06d}"          # Archive naming format
    eng.io["output_freq"] = 1000                 # Jobs per archive file

**Parameters:**

- **output**: String specifying the base directory name for storing output files
- **output_format**: Format string for archive file naming (supports Python string formatting)
- **output_freq**: Integer defining how many completed jobs to include per archive file

External Program Paths
~~~~~~~~~~~~~~~~~~~~~~

Specify paths to required external software:

.. code-block:: python

    # External program configuration
    eng.programs["mgltools"] = "/path/to/MGLTools-1.5.6"           # MGLTools installation
    eng.programs["vina"] = "/path/to/vina"                         # AutoDock Vina executable
    eng.programs["vinardo"] = "/path/to/smina"                     # Smina with Vinardo scoring (optional)
    eng.programs["voromqa"] = "/path/to/voronota-voromqa"          # VoroMQA executable

**Required Programs:**

- **mgltools**: Path to MGLTools installation directory (required for PDBQT conversion)
- **vina**: Path to AutoDock Vina executable (required for molecular docking)
- **vinardo**: Path to Smina executable with Vinardo scoring (optional alternative scorer)
- **voromqa**: Path to VoroMQA executable (optional for advanced model quality assessment)

Variant Generation Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Control how protein variants are generated:

.. code-block:: python

    # Variant generation settings
    eng.variant["name"] = "mutation"                   # Generation strategy
    eng.variant["selection"] = "A:100 A:150 A:200"     # Residues to mutate
    eng.variant["fixed"] = "A:50 A:75"                 # Fixed residues during optimization
    eng.variant["excluded_all"] = "CGP"                # Globally excluded amino acids
    eng.variant["excluded"] = {"A:100": "FWYM"}        # Position-specific exclusions
    eng.variant["conservative"] = True                 # Use conservative mutations
    eng.variant["max_iterations"] = 1000               # Maximum variants to generate
    eng.variant["combinations"] = 2                    # Maximum simultaneous mutations
    eng.variant["matrix"] = "blosum62"                 # Substitution matrix
    eng.variant["msa"] = "sequences.fasta"             # Multiple sequence alignment file

**Strategy Options:**

- **name**:
  - ``"mutation"``: Matrix-based mutagenesis using substitution matrices
  - ``"exhaustive"``: Systematic combinatorial mutagenesis
  - ``"msa"``: Sequence variants from FASTA files

**Selection Parameters:**

- **selection**: Space-separated list of residues to mutate (format: "ChainID:ResidueNumber")
- **fixed**: Space-separated list of residues to keep fixed during MODELLER optimization (format: "ChainID:ResidueNumber")
- **excluded_all**: String of amino acid single-letter codes to exclude globally
- **excluded**: Dictionary mapping residue positions to excluded amino acids

**Mutation Control:**

- **conservative**: Boolean controlling mutation bias (True = conservative, False = non-conservative)
- **max_iterations**: Maximum number of variants to generate
- **combinations**: Maximum number of simultaneous mutations per variant
- **matrix**: Substitution matrix name ("blosum62" or custom matrix specification)
- **msa**: Path to FASTA file containing sequences for FASTA-based variant generation

Structure Modeling Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure 3D structure modeling with MODELLER:

.. code-block:: python

    # MODELLER configuration
    eng.model["name"] = "modeller"             # Modeling method
    eng.model["num_models"] = 5                # Number of models per variant
    eng.model["optimize_radius"] = 8           # Optimization radius in Angstroms
    eng.model["optimize_level"] = 1            # Optimization thoroughness

**Parameters:**

- **name**: Modeling method (currently only "modeller" is supported)
- **num_models**: Number of 3D models to generate per sequence variant
- **optimize_radius**: Distance in Angstroms around mutations to optimize (0 = just mutated residues)
- **optimize_level**:
  - ``0``: Fast optimization (very_fast schedule, fast MD)
  - ``1``: Normal optimization (normal schedule, slow MD)
  - ``2``: Thorough optimization (slow schedule, very_slow MD)

Model Quality Assessment Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure quality filtering for generated models:

.. code-block:: python

    # Model quality thresholds
    eng.model_quality["norm_dope"] = -1.0      # Normalized DOPE score threshold
    eng.model_quality["voromqa"] = 0.4         # VoroMQA score threshold

**Quality Metrics:**

- **norm_dope**: Normalized DOPE score threshold (models with scores higher than this value are rejected)
- **voromqa**: VoroMQA score threshold (models with scores below this value are rejected, requires VoroMQA installation)

Molecular Docking Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure protein-ligand docking parameters:

.. code-block:: python

    # Docking configuration
    eng.evaluator["name"] = "vina"                         # Docking engine
    eng.evaluator["exhaustiveness"] = 32                   # Search thoroughness
    eng.evaluator["box_center"] = [10.0, 15.0, 20.0]       # Search box center (x, y, z)
    eng.evaluator["box_size"] = [20.0, 20.0, 20.0]         # Search box dimensions

**Docking Parameters:**

- **name**:
  - ``"vina"``: Standard AutoDock Vina scoring
  - ``"vinardo"``: Vinardo scoring function (using Smina)
- **exhaustiveness**: Search thoroughness (higher values = more thorough but slower)
- **box_center**: List of three floats defining the center coordinates of the docking search box
- **box_size**: List of three floats defining the dimensions of the search box in Angstroms

Ligand and Measurement Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure ligands and distance measurements:

.. code-block:: python

    # Add ligands and define measurements
    ligand = eng.add_ligand("substrate", "ligand.pdbqt")

    # Distance measurements between protein and ligand atoms
    ligand.add_measurement("catalytic_distance", "distance",
                          "chainID A and resid 87 and name N", "name O10")

**Ligand Methods:**

- **add_ligand(name, filename)**: Add a ligand for docking
  - ``name``: String identifier for the ligand
  - ``filename``: Path to PDBQT file containing ligand structure

**Measurement Methods:**

- **add_measurement(name, metric, protein_selection, ligand_selection)**: Add distance measurement
  - ``name``: Unique identifier for the measurement
  - ``metric``: Currently only "distance" is supported
  - ``protein_selection``: MDAnalysis selection string for protein atoms
  - ``ligand_selection``: MDAnalysis selection string for ligand atoms

Selection Syntax Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~

GDEE uses MDAnalysis selection syntax for specifying atoms. Common examples:

.. code-block:: python

    # Atom-based selections
    "chainID A and resid 100 and name CA"     # Atom CA of residue 100 in chain A
    "name O10"                                # Atom O10 of the ligand

For complete syntax documentation, refer to the MDAnalysis selection documentation.


Example Workflows
-----------------

Gene Discovery
~~~~~~~~~~~~~~

Gene discovery workflows focus on exploring sequence variants from FASTA files to identify naturally occurring variants that may catalyze specific biochemical reactions.

.. code-block:: python

    from gdee import ProteinEngineering

    # Initialize the GDEE platform with protein name and database file
    # This creates the main workflow orchestrator and database connection
    eng = ProteinEngineering("target-protein", "db-filename")

    # Set the template PDB structure for homology modeling
    # This structure serves as the basis for generating 3D models of variants
    eng.pdb = "protein.pdb"

    # Configure paths to external programs required for the workflow
    # MGLTools: Used for converting PDB to PDBQT format for docking
    eng.programs["mgltools"] = "/path/to/MGLTools-1.5.6"
    # AutoDock Vina: Molecular docking engine for protein-ligand interactions
    eng.programs["vina"] = "/path/to/autodock_vina/vina"
    # VoroMQA: Model quality assessment using Voronoi tessellation
    eng.programs["voromqa"] = "/path/to/voronota/voronota-voromqa"

    # Configure file archiving and output management
    # Archive successful jobs to compressed files for storage efficiency
    eng.io["output"] = "files"
    # Naming format for archived files (5-digit zero-padded numbers)
    eng.io["output_format"] = ".{:05d}"
    # Archive every 50 completed jobs to balance I/O and storage
    eng.io["output_freq"] = 50

    # Configure distributed execution platform
    # Use MPI for parallel processing across multiple compute nodes
    eng.platform["name"] = "mpi"
    # Number of CPU cores to use on each MPI node for local parallelism
    eng.platform["local_cpu"] = 16

    # Configure 3D structure modeling parameters
    # Number of models to be used in the docking step
    eng.model["num_models"] = 5
    # Optimize residues within 8 Å of mutation sites (local optimization)
    eng.model["optimize_radius"] = 8
    # Use moderate optimization level (0=fast, 1=normal, 2=slow)
    eng.model["optimize_level"] = 1

    # Configure molecular docking parameters
    # Use AutoDock Vina as the docking engine
    eng.evaluator["name"] = "vina"
    # High exhaustiveness for thorough conformational search (default: 8)
    eng.evaluator["exhaustiveness"] = 200
    # Center coordinates of the docking search box (x, y, z in Angstroms)
    eng.evaluator["box_center"] = [32, 20, 39.5]
    # Dimensions of the search box (width, height, depth in Angstroms)
    eng.evaluator["box_size"] = [14, 11, 16]

    # Add a ligand for docking
    # Returns a Ligand object for further configuration
    lig = eng.add_ligand("ligand-name", "ligand-file.pdbqt")

    # Define distance measurements between protein and ligand atoms
    # These measurements will be computed for each docking pose
    # Measure distance from residue 292 CA atom to ligand CL1 atom
    lig.add_measurement("metric-1", "distance", "resid 292 and name CA", "name CL1")
    # Measure distance from residue 160 N atom to ligand O2 atom
    lig.add_measurement("metric-2", "distance", "resid 160 and name N", "name O2")
    # Measure distance from residue 220 NH2 group to ligand C12 atom
    lig.add_measurement("metric-3", "distance", "resid 220 and name NH2", "name C12")

    # Configure variant generation from FASTA file
    print("Running MSA")
    eng.variant["name"] = "msa"
    # FASTA file containing sequences to be evaluated
    eng.variant["msa"] = "sequences_20_identity.fasta"

    # Execute the complete gene discovery workflow
    # This will: generate variants → model structures → assess quality → dock ligands → measure distances
    eng.run()

    print("Done!")



Enzyme Engineering
~~~~~~~~~~~~~~~~~~

Enzyme engineering workflows focus on mutation-based optimization to improve specific biochemical reactions.

.. code-block:: python

    from gdee import ProteinEngineering

    # Initialize platform for enzyme engineering study
    eng = ProteinEngineering("target-protein", "db-filename")
    # Template structure of the enzyme to be engineered
    eng.pdb = "protein.pdb"

    # Configure external program paths (same as gene discovery)
    eng.programs["mgltools"] = "/path/to/MGLTools-1.5.6"
    eng.programs["vina"] = "/path/to/autodock_vina/vina"
    eng.programs["voromqa"] = "/path/to/voronota/voronota-voromqa"

    # Configure archiving for larger datasets typical in enzyme engineering
    eng.io["output"] = "files"
    eng.io["output_format"] = ".{:05d}"
    # Archive less frequently due to larger number of variants
    eng.io["output_freq"] = 2000

    # Use distributed execution for systematic mutagenesis screens
    eng.platform["name"] = "mpi"
    eng.platform["local_cpu"] = 16

    # Define structural constraints for mutagenesis
    # Fixed positions that will not be optimized by MODELLER (critical for structure/function)
    # Format: "ChainID:ResidueNumber" separated by spaces
    eng.variant["fixed"] = "A:237 A:206 A:160 A:87 A:161"
    # Positions selected for mutagenesis (target sites for optimization)
    eng.variant["selection"] = "A:185 A:159 A:88 A:93 A:92 A:208"

    # Define amino acid groups for systematic exclusion rules
    # Special amino acids that might disrupt structure
    SPECIAL = "CGP"  # Cysteine (disulfide), Glycine (flexible), Proline (rigid)
    # Charged amino acids
    POSITIVE = "RHK"  # Arginine, Histidine, Lysine
    NEGATIVE = "DE"   # Aspartate, Glutamate
    # Polar amino acids
    POLAR = "QNST"    # Glutamine, Asparagine, Serine, Threonine
    # Hydrophobic amino acids
    NON_POLAR = "AVILMFYW"  # Ala, Val, Ile, Leu, Met, Phe, Trp, Tyr

    # Global exclusion rule: don't use these amino acids at any position
    eng.variant["excluded_all"] = NEGATIVE + POSITIVE + SPECIAL

    # Position-specific exclusion rules for fine-tuned mutagenesis
    # Each position has customized restrictions based on structural role
    eng.variant["excluded"] = {
        "A:159": POLAR,      # Exclude polar residues at position 159
        "A:88": "FWY",       # Exclude bulky aromatics at position 88
        "A:93": "FWYM",      # Exclude aromatics and Met at position 93
        "A:92": NON_POLAR,   # Exclude hydrophobic residues at position 92
        "A:208": "FWYM",     # Exclude bulky residues at position 208
    }

    # Number of models to be used in the docking step
    eng.model["num_models"] = 5
    # Optimize residues within 8 Å of mutation sites (local optimization)
    eng.model["optimize_radius"] = 8
    # Use moderate optimization level (0=fast, 1=normal, 2=slow)
    eng.model["optimize_level"] = 1

    # Configure molecular docking parameters
    # Name of the docking engine to use
    eng.evaluator["name"] = "vina"
    # Exhaustiveness value used in the docking calculations
    eng.evaluator["exhaustiveness"] = 200
    # Center coordinates of the docking search box (x, y, z in Angstroms)
    eng.evaluator["box_center"] = [32, 20, 39.5]
    # Dimensions of the search box (width, height, depth in Angstroms)
    eng.evaluator["box_size"] = [14, 11, 16]

    # Add a ligand for docking
    # Returns a Ligand object for further configuration
    lig = eng.add_ligand("ligand-name", "ligand-file.pdbqt")

    # Define critical distance measurements for pose filtering
    # Distance from catalytic residue N87 to substrate O10 (catalytic interaction)
    lig.add_measurement("metric1", "distance", "chainID A and resid 87 and name N", "name O10")
    # Distance from binding residue N161 to substrate O10 (substrate positioning)
    lig.add_measurement("metric2", "distance", "chainID A and resid 161 and name N", "name O10")
    # Distance from S160 side chain to substrate C3 (cofactor interaction)
    lig.add_measurement("metric3", "distance", "chainID A and resid 160 and name OG", "name C3")
    lig.add_measurement("metric4", "distance", "chainID A and resid 185 and name CA", "name C2")

For enzyme engineering, two variant generation strategies are available:

1. Exhaustive combinatorial mutagenesis: Screen all possible double amino acid substitutions in selected residues

.. code-block:: python

    # Exhaustive combinatorial mutagenesis
    # Generate all possible combinations at selected positions
    print("Running all 2 by 2")
    eng.variant["name"] = "exhaustive"
    # Maximum number of simultaneous mutations per variant (In this example will generate single and double mutamts)
    eng.variant["combinations"] = 2

    # Execute pipeline
    eng.run()

    print("Done!")

2. Sampling-based mutagenesis: Sample the mutation space using conservative substitutions and scoring matrices

.. code-block:: python

    # Sampling approach using substitution matrices (e.g., BLOSUM62) to prioritize conservative changes
    print("Sampling 2 by 2")
    eng.variant["name"] = "mutation"
    # Use conservative amino acid substitutions based on evolutionary data
    eng.variant["conservative"] = True
    # Maximum number of variants to generate (prevents excessive library sizes)
    eng.variant["max_iterations"] = 50000
    # Maximum number of simultaneous mutations per variant (In this example will generate single and double mutamts)
    eng.variant["combinations"] = 2

    # Execute pipeline
    eng.run()

    print("Done!")


Filtering Docking Poses and Ranking Variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After workflow completion, use the analysis module to identify the most promising variants based on multiple criteria:

.. code-block:: python

    from gdee.database import Database
    from gdee.analysis.filters import Metric

    # Connect to the results database generated by the platform
    db = Database("db-filename.sqlite3")

    # Create metric objects for each measurement defined in the workflow
    # These correspond to the distance measurements added to ligands above
    metric1 = Metric("metric1", db)
    metric2 = Metric("metric2", db)
    metric3 = Metric("metric3", db)
    metric4 = Metric("metric4", db)

    # Define filtering rules based on biochemical criteria
    # Close catalytic contact for example (< 3.5 Å indicates proper positioning)
    rule1 = metric1 < 3.5
    # Close binding interaction for example (< 3.5 Å indicates strong binding)
    rule2 = metric2 < 3.5
    # Cofactor should be closer to substrate than pocket residue (selectivity)
    rule3 = metric4 > metric3

    # Combine rules using boolean logic
    # Variants must have either good catalytic OR binding contacts AND proper selectivity
    rule4 = (rule1 | rule2) & rule3

    # Create ranking based on binding energy (lower is better)
    # Variants passing the filtering rules are ranked by docking energy
    rank = rule4.rank()
    # Export all passing variants to SQLite database for further analysis
    rank.export_sqlite("rank.sqlite3", "All_Variants")

    # Filter by mutation count for systematic analysis
    # Single mutants are often easier to validate experimentally
    single = rank.by_num_mutations(1)
    single.export_sqlite("rank.sqlite3", "Single")
    # Export top 15 single mutants to CSV for easy inspection
    single.export_csv("rank_single.csv", 15)

    # Double mutants may show synergistic effects
    double = rank.by_num_mutations(2)
    double.export_sqlite("rank.sqlite3", "Double")
    # Export top 15 double mutants for experimental follow-up
    double.export_csv("rank_double.csv", 15)


Rescoring with Machine Learning Scoring Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After completing initial docking and ranking workflows, you can rescore the candidates using a machine learning-based scoring function.

.. note::

   Machine learning scoring functions must be trained prior to use. The rescoring workflow requires pre-trained model files in pickle format:

   - **RF-Score models**: Random Forest-based scoring functions (v1, v2, v3 variants)
   - **PLECnn model**: Deep learning model using Protein-Ligand Extended Connectivity interaction fingerprints

For training these models, refer to the ODDT package documentation.


.. code-block:: python

    from gdee.database import Ranked_Database
    from gdee.engineer import RescoreVariants

    # Load the ranked database from previous workflow filtering
    input_db = Ranked_Database("rank.sqlite3")

    # Specify output database for rescored results
    # This will create a new database with all rescored poses
    output_db = "rank_rescored.sqlite3"

    # Initialize rescoring workflow
    rescoring = RescoreVariants(
        input_db,                                       # Input database with ranked variants
        output_db,                                      # Output database for rescored results
        "AllVariants",                                  # Table name containing variants to rescore
        [
            "RFScore_v1_pdbbind2016.pickle",            # path to RF-Score v1 model
            "RFScore_v2_pdbbind2016.pickle",            # path to RF-Score v2 model
            "RFScore_v3_pdbbind2016.pickle",            # path to RF-Score v3 model
            "PLECnn_p5_l1_pdbbind2016_s65536.pickle"    # path to PLECnn deep learning model
        ],
        "files"                                         # Directory containing structure files
    )

    # Configure distributed or simple execution for rescoring
    rescoring.platform["name"] = "mpi"
    # Number of CPU cores to use on each MPI node for local parallelism
    rescoring.platform["local_cpu"] = 16

    # Execute rescoring workflow
    print("Running rescore")
    rescoring.run()
    print("Done!")



Useful Scripts
--------------

Script to run BLAST search against Swissprot database and save results in XML format:

.. code-block:: python

    from Bio import SeqIO
    from Bio.Blast import NCBIWWW

    # Configuration parameters for BLAST search
    FASTA = "target-protein.fasta"        # Input protein sequence file
    DATABASE = "swissprot"                # database
    E_VALUE = 10                          # E-value threshold
    MAX_HITS = 10000                      # Maximum number of sequences to retrieve
    MIN_IDENTITY = 10                     # Minimum percent identity threshold
    OUTPUT = "blast_results.xml"          # Output file in XML format for parsing

    # Read the query protein sequence from FASTA file
    query = SeqIO.read(FASTA, format="fasta")

    # Perform online BLAST search
    result = NCBIWWW.qblast("blastp", DATABASE, query.seq, expect=E_VALUE,
                            hitlist_size=MAX_HITS, perc_ident=MIN_IDENTITY)

    # Save BLAST results in XML format for subsequent parsing
    # XML format preserves all alignment details and statistics
    with open(OUTPUT, "w") as fd:
        fd.write(result.read())


Script to filter BLAST results based on coverage and identity, and save sequences in FASTA format:

.. code-block:: python

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Blast import NCBIXML

    # Configuration parameters for filtering BLAST results
    OUTPUT = "sequences_single_20_identity.fasta"   # Output FASTA file with filtered sequences
    BLAST_XML = "blast_results.xml"                 # Input XML file from BLAST search
    MIN_COVERAGE = 0.8                              # Minimum query coverage
    MIN_IDENTITY = 0.2                              # Minimum sequence identity

    # Parse BLAST results from XML file
    with open("blast_results.xml") as fd:
        blast = NCBIXML.read(fd)

    # Initialize lists to collect filtered sequences and statistics
    expects = []      # E-values for statistical analysis
    sequences = []    # Filtered sequences for output

    # Iterate through all database matches (alignments) from BLAST
    for aln in blast.alignments:
        # Get the best HSP (High-scoring Segment Pair) for each alignment
        # HSPs represent local alignments between query and database sequence
        hsp = aln.hsps[0]

        # Calculate query coverage as fraction of query sequence aligned
        # Higher coverage indicates more complete homology
        coverage = (hsp.query_end - hsp.query_start + 1) / blast.query_length

        # Calculate sequence identity as fraction of identical residues
        # Higher identity indicates closer evolutionary relationship
        identity = hsp.identities / hsp.align_length

        # Skip identical sequences (identity = 1.0) to avoid self-matches
        # These don't provide new information for diversity analysis
        if identity == 1:
            continue

        # Apply coverage and identity filters to select high-quality homologs
        # These thresholds ensure sequences are sufficiently similar and complete
        if coverage >= MIN_COVERAGE and identity >= MIN_IDENTITY:
            # Extract aligned subject sequence and remove gap characters
            # Gaps (-) are alignment artifacts and should be removed
            seq = Seq(hsp.sbjct).replace("-", "")

            # Create sequence record with database accession as identifier
            # Empty description and annotation fields for simplicity
            record = SeqRecord(seq, aln.accession, "", "")
            sequences.append(record)

            # Store E-value for statistical summary
            expects.append(hsp.expect)

    # Print summary statistics about filtering results
    # Helps assess the quality and diversity of the filtered dataset
    print("Filtered {} sequences from {}".format(len(sequences), len(blast.alignments)))
    print("Expect values. Max: {}. Min: {}".format(max(expects), min(expects)))

    # Write filtered sequences to FASTA file
    SeqIO.write(sequences, OUTPUT, "fasta")


The resulting FASTA file is ready for use in the GDEE platform for gene discovery workflows.
