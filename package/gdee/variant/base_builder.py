"""
Base class for protein variant builders in the GDEE platform.

This module provides the abstract base class that defines the common interface
and functionality for all variant generation strategies.
"""


from .sequence import ProtSeq


class BaseBuilder:
    """
    Abstract base class for protein variant generators.
    
    This class defines the common interface and shared functionality for all
    variant builders including database interaction, exclusion handling,
    and result management. Concrete implementations include MSABuilder,
    MutationBuilder, and ExhaustiveBuilder.
    
    Attributes:
        parameters (dict): Configuration parameters from ProteinEngineering
        db (Database): SQLite database connection for result storage
        prot_id (int): Database protein ID
        protein (ProtSeq): Wildtype protein sequence
        _initialized (bool): Lazy initialization flag
        _excluded (dict): Per-residue amino acid exclusion rules
        _variants (set): Cache of existing variant names
    """
    
    def __init__(self, parameters, database):
        """
        Initialize base builder with configuration and database.
        
        Args:
            parameters (dict): Builder configuration parameters
            database (Database): Database connection for result storage
        """
        self.parameters = parameters
        self.db = database
        self._initialized = False
        self.prot_id = None
        self.protein = None
        self._excluded = parameters["excluded"]
        self._variants = set()

    def is_excluded(self, residue, code):
        """
        Check if amino acid is excluded for a specific residue.
        
        Args:
            residue (SeqPos): Target residue position
            code (str): Single-letter amino acid code
            
        Returns:
            bool: True if amino acid is excluded for this position
        """
        key = "{}:{}".format(residue.chain, residue.resid)
        return code in self._excluded.get(key, "")

    def initialize(self):
        """
        Perform lazy initialization of builder components.
        
        Registers protein in database, loads sequence from PDB file,
        caches existing variants, and calls strategy-specific initialization.
        """
        self._initialized = True
        self.prot_id = self.db.register_protein(self.parameters["protein_name"])
        self.protein = ProtSeq(self.parameters["protein_name"], self.parameters["pdb_file"])
        # Faster than making queries and low memory overhead
        self._variants.update(item[0] for item in self.db.fetch_variants(self.prot_id))
        self.special_initialize()

    def special_initialize(self):
        """
        Strategy-specific initialization (abstract method).
        
        Must be implemented by concrete builder classes to set up
        strategy-specific data structures and parameters.
        
        Raises:
            NotImplementedError: If not implemented by child class
        """
        raise NotImplementedError("Child classes must implement this method")

    def next_job(self):
        """
        Get the next variant generation job.
        
        Handles lazy initialization and variant name registration.
        This is the main entry point called by Pipeline.next_job().
        
        Returns:
            DataContainer or None: Job data for pipeline processing,
                                   or None if no more variants available
        """
        if not self._initialized:
            self.initialize()  # Lazy initialization

        job = self.fetch_next_job()
        if job is not None:
            self.add_variant(job.variant.name)

        return job

    def fetch_next_job(self):
        """
        Generate the next variant (abstract method).
        
        Must be implemented by concrete builder classes to provide
        strategy-specific variant generation logic.
        
        Returns:
            DataContainer or None: Job data container or None if finished
            
        Raises:
            NotImplementedError: If not implemented by child class
        """
        raise NotImplementedError("Child classes must implement this method")

    def variant_exists(self, name):
        """
        Check if a variant name already exists.
        
        Args:
            name (str): Variant name to check
            
        Returns:
            bool: True if variant already exists
        """
        return name in self._variants

    def add_variant(self, name):
        """
        Register a new variant name.
        
        Args:
            name (str): Variant name to register
            
        Raises:
            RuntimeError: If variant name already exists
        """
        if name in self._variants:
            raise RuntimeError("Variant {} already exists".format(name))

        self._variants.add(name)

    def unsave_results(self, data):
        """
        Remove variant results from database (error recovery).
        
        Called when pipeline processing fails to clean up partial results.
        
        Args:
            data (DataContainer): Job data with variant information
        """
        name = data.variant.name
        self.db.remove_variant(name)
        self._variants.discard(name)

    def save_results(self, data):
        """
        Save variant processing results to database.
        
        This method orchestrates the storage of all pipeline results including
        variant information, 3D models, quality scores, docking evaluations,
        poses, and measurements.
        
        Args:
            data (DataContainer): Complete job results from pipeline processing
            
        Returns:
            bool: True if results were saved successfully
        """
        name = data.variant.name

        # Sanity check for processing completeness
        has_eval = False
        all_rejected = True
        try:
            for model in data.modeling.models:
                if model.evals:
                    has_eval = True

                if not model.rejected:
                    all_rejected = False
        except AttributeError:
            pass

        if data.fatal_error:
            print("Error while processing variant: {}".format(name))

        if all_rejected:
            print("Warning: all models rejected during quality assessment for variant '{}'.".format(name))

        if not has_eval:
            print("Error: no evaluations for variant: {}".format(name))

        if data.fatal_error or not has_eval or all_rejected:
            return False

        # Register variant in database
        variant_id = self.db.register_variant(self.prot_id, name,
                                              data.variant.to_modeller(),
                                              data.variant_dir,
                                              data.is_wildtype)

        # Save modeling results
        modeling = data.modeling
        for model in modeling.models:
            if not model.evals:
                continue

            # Register 3D model with quality scores
            model_id = self.db.register_model(variant_id, modeling.method,
                                              model.scores.jsonfy(), model.pdb,
                                              model.rejected)

            # Save docking evaluations for each ligand
            for ligand_name, evaluation in model.evals.items():
                eval_id = self.db.register_evaluation(variant_id, model_id, evaluation)
                pose_id_list = self.db.register_poses(eval_id, evaluation.energies)

                # Save measurement data for pose analysis
                self.db.register_measurements(eval_id, pose_id_list,
                                              evaluation.measurements)

        print("Ended with variant '{}'".format(data.variant.name))
        return True
