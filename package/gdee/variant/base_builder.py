"""Base class for protein variant generation."""


from .sequence import ProtSeq


class BaseBuilder:
    """Base class for variant generation strategies.

    Handles database interaction, exclusion rules, and result storage.
    Subclasses implement specific variant generation strategies.
    """

    def __init__(self, parameters, database):
        """Initialize variant builder.

        Args:
            parameters: Configuration dictionary with builder parameters
            database: Database connection instance
        """
        self.parameters = parameters
        self.db = database
        self._initialized = False
        self.prot_id = None
        self.protein = None
        self._excluded = parameters["excluded"]
        self._variants = set()

    def is_excluded(self, residue, code):
        """Check if amino acid is excluded at a position.

        Args:
            residue: SeqPos object representing position
            code: Single-letter amino acid code

        Returns:
            bool: True if amino acid is excluded
        """
        key = "{}:{}".format(residue.chain, residue.resid)
        return code in self._excluded.get(key, "")


    def initialize(self):
        """Perform lazy initialization of protein and database.

        Loads protein sequence and registers in database.
        """
        self._initialized = True
        self.prot_id = self.db.register_protein(self.parameters["protein_name"])
        self.protein = ProtSeq(self.parameters["protein_name"],
                               self.parameters["pdb_file"])
        # Faster than making queries and low memory overhead
        self._variants.update(item[0] for item in self.db.fetch_variants(self.prot_id))
        self.special_initialize()

    def special_initialize(self):
        """Subclass-specific initialization.

        Must be implemented by subclasses.

        Raises:
            NotImplementedError: Must be implemented by subclass
        """
        raise NotImplementedError("Child classes must implement this method")

    def next_job(self):
        """Get next variant generation job.

        Returns:
            DataContainer: Job data with variant information, or None if complete
        """
        if not self._initialized:
            self.initialize()

        job = self.fetch_next_job()
        if job is not None:
            self.add_variant(job.variant.name)

        return job

    def fetch_next_job(self):
        """Generate next variant job.

        Must be implemented by subclasses.

        Returns:
            DataContainer: Job data or None if generation complete

        Raises:
            NotImplementedError: Must be implemented by subclass
        """
        raise NotImplementedError("Child classes must implement this method")

    def variant_exists(self, name):
        """Check if variant already exists.

        Args:
            name: Variant name

        Returns:
            bool: True if variant exists
        """
        return name in self._variants

    def add_variant(self, name):
        """Register a generated variant.

        Args:
            name: Variant name

        Raises:
            RuntimeError: If variant already exists
        """
        if name in self._variants:
            raise RuntimeError("Variant {} already exists".format(name))
        self._variants.add(name)

    def unsave_results(self, data):
        """Remove variant from database on failure.

        Args:
            data: Job data with variant information
        """
        name = data.variant.name
        variant_id = self.db.remove_variant(self.prot_id, name)
        self._variants.add(name)

    def save_results(self, data):
        """Save variant processing results to database.

        Stores variant sequence, models, docking results, and measurements.

        Args:
            data: Job results with variant, models, and evaluations

        Returns:
            bool: True if saved successfully, False if validation failed
        """
        name = data.variant.name

        # Sanity check
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
            return

        variant_id = self.db.register_variant(self.prot_id, name,
                                              data.variant.to_modeller(),
                                              data.variant_dir,
                                              data.is_wildtype)

        modeling = data.modeling
        for model in modeling.models:
            if not model.evals:
                continue

            model_id = self.db.register_model(variant_id, modeling.method,
                                              model.scores.jsonfy(), model.pdb,
                                              model.rejected)

            for ligand_name, evaluation in model.evals.items():
                eval_id = self.db.register_evaluation(variant_id, model_id,
                                                      evaluation)
                pose_id_list = self.db.register_poses(eval_id, evaluation.energies)

                self.db.register_measurements(eval_id, pose_id_list,
                                              evaluation.measurements)

        print("Ended with variant '{}'".format(data.variant.name))
        return True
