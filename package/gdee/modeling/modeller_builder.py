"""
MODELLER-based 3D structure modeling for protein variants.

This module provides homology modeling functionality using MODELLER to generate
3D structures for protein variants with mutation-specific optimization strategies.
"""

from path import Path
import modeller as mdl
import MDAnalysis as mda
from modeller import automodel
import shutil
import contextlib
from tempfile import TemporaryDirectory
from gdee.misc import DataContainer
import warnings

warnings.filterwarnings("ignore", module=r"MDAnalysis.*")


class ModellerBuilder:
    """
    MODELLER-based 3D structure builder for protein variants.
    
    This class generates 3D structural models for protein variants using
    comparative modeling. It supports mutation-specific optimization with
    configurable optimization levels and local refinement strategies.
    
    Attributes:
        parameters (dict): Modeling configuration including template PDB,
                          optimization settings, and model generation parameters
        env (modeller.environ): MODELLER environment for model building
    """
    
    def __init__(self, parameters):
        """
        Initialize MODELLER environment and configuration.
        
        Sets up MODELLER environment with optimized parameters for
        GDEE workflows including HETATM and water handling.
        
        Args:
            parameters (dict): Configuration including template PDB file,
                             optimization level, and model count
        """
        self.parameters = parameters
        
        # Configure MODELLER environment
        mdl.log.none()  # Suppress verbose output
        self.env = mdl.environ()
        self.env.io.hetatm = True      # Include heteroatoms (ligands, cofactors)
        self.env.io.water = True       # Include water molecules
        self.env.edat.dynamic_lennard = True  # Dynamic Lennard-Jones potential
        self.env.schedule_scale[mdl.physical.lennard_jones] = 1.0
        mdl.log.level(0, 0, 0, 0, 0)   # Minimal logging

    def run(self, job_data):
        """
        Generate 3D models for a protein variant.
        
        This method:
        1. Creates temporary workspace with template structure
        2. Writes alignment file for variant vs template
        3. Builds multiple models using MODELLER
        4. Renames residues to match variant sequence
        5. Saves models to job directory

        Args:
            job_data (DataContainer): Job data with variant information
            
        Returns:
            DataContainer: Updated job data with modeling results
        """
        temp_dir = TemporaryDirectory(prefix="gdee_modeller")
        temp_path = Path(temp_dir.name)
        Path(self.parameters["pdb_file"]).copy(temp_path / "template.pdb")

        with temp_path:
            # Create alignment file for comparative modeling
            self.write_alignment(job_data)

            # Build multiple models with MODELLER
            model_data = self.build_models(job_data)
            if not model_data:
                job_data.fatal_error = True
                return job_data

            # Select best models based on Normalized DOPE scores
            model_data.sort(key=lambda x: x["Normalized DOPE score"])
            top_models = model_data[:self.parameters["num_models"]]

            # Load models into MDAnalysis for processing
            pdb_list = [str(model["name"]) for model in top_models]
            structure = mda.Universe(pdb_list[0], pdb_list)
            self.rename_models(structure, job_data.variant)

            # Save processed models to job directory
            model_list = []
            for i, ts in enumerate(structure.trajectory):
                filename = "model_{:04d}.pdb".format(i)
                structure.atoms.write(str(job_data.job_dir / filename))

                # Create model data container
                model = DataContainer()
                model.scores = DataContainer()
                model.scores.norm_dope = top_models[i]["Normalized DOPE score"]
                model.scores.molpdf = top_models[i]["molpdf"]
                model.pdb = filename
                model.rejected = False  # Will be set by quality checker
                model_list.append(model)

        # Store modeling results in job data
        if model_list:
            modeled = DataContainer()
            modeled.method = "modeller"
            modeled.models = model_list
            job_data.modeling = modeled
        else:
            job_data.fatal_error = True

        return job_data

    def write_alignment(self, job_data):
        """
        Create PIR alignment file for MODELLER.
        
        Generates alignment between template structure and variant sequence
        in PIR format required by MODELLER for comparative modeling.
        
        Args:
            job_data (DataContainer): Job data with wildtype and variant sequences
        """
        template = (">P1;template\n"
                   "structureX:template.pdb:.:.:.:.::::\n"
                   "{}*\n\n"
                   ">P1;model\n"
                   "sequence:model.pdb:.:.:.:.::::\n"
                   "{}*\n")
        
        with open("alignment.ali", "w") as fd:
            fd.write(template.format(
                job_data.wildtype.to_modeller(),
                job_data.variant.to_modeller()
            ))

    def build_models(self, job_data):
        """
        Build 3D models using MODELLER with mutation-specific optimization.
        
        Configures MODELLER parameters based on optimization level and
        builds multiple models with assessment scores.
        
        Args:
            job_data (DataContainer): Job data with mutation information
            
        Returns:
            list: List of model data dictionaries with scores and filenames
        """
        # Create MODELLER automodel instance
        model = MutationModel(
            self.env,
            alnfile="alignment.ali",
            knowns="template",
            sequence="model",
            assess_methods=(automodel.assess.normalized_dope),
        )

        # Generate multiple models for selection
        model.starting_model = 1
        model.ending_model = 3 * self.parameters["num_models"]  # Build extra for selection
        model.final_malign3d = True  # Final alignment optimization

        # Configure optimization level
        opt_level = self.parameters["optimize_level"]
        if opt_level == 0:
            # Fast optimization for development/testing
            model.library_schedule = automodel.autosched.very_fast
            model.md_level = automodel.refine.fast
            model.max_var_iterations = 100
            model.repeat_optimization = 1

        elif opt_level == 1:
            # Normal optimization for standard use
            model.library_schedule = automodel.autosched.normal
            model.md_level = automodel.refine.slow
            model.max_var_iterations = 200
            model.repeat_optimization = 1

        else:
            # Slow optimization for high-accuracy modeling
            model.library_schedule = automodel.autosched.slow
            model.md_level = automodel.refine.very_slow
            model.max_var_iterations = 300
            model.repeat_optimization = 1

        # Build models with mutation-specific optimization
        with contextlib.redirect_stdout(None):  # Suppress MODELLER output
            model.select_opt_residues(
                job_data.mut_index,      # Mutated residue indices
                job_data.fixed_index,    # Fixed residue indices
                self.parameters["optimize_radius"]  # Optimization sphere radius
            )
            
            model.make()  # Execute model building

        # Extract model data with scores
        model_data = []
        for data in model.outputs:
            if "Normalized DOPE score" not in data:
                continue  # Skip models without quality scores

            data["name"] = Path(data["name"]).stem + "_fit.pdb"
            model_data.append(data)

        return model_data

    def rename_models(self, structure, prot_seq):
        """
        Rename residues in models to match variant sequence.
        
        Updates residue names, numbers, and chain assignments to match
        the variant sequence for consistency with downstream analysis.
        
        Args:
            structure (mda.Universe): MDAnalysis universe with model structures
            prot_seq (ProtSeq): Variant protein sequence for renaming
        """
        res_idx = 0
        for chain in prot_seq:
            # Create new segment for each chain
            segment = structure.add_Segment(segid=chain.code)

            for seq_pos in chain:
                if seq_pos.is_gap:
                    continue  # Skip gaps in sequence

                # Update residue properties
                residue = structure.residues[res_idx]
                residue.resid = seq_pos.resid
                residue.resname = seq_pos.resname
                residue.segment = segment
                res_idx += 1


class MutationModel(automodel.automodel):
    """
    Custom MODELLER automodel class with mutation-specific optimization.
    
    Extends MODELLER's automodel to support selective optimization of
    mutated and nearby residues while keeping distant regions fixed.
    This improves efficiency and focuses refinement on relevant areas.
    """
    
    def select_opt_residues(self, residues, excluded, radius):
        """
        Configure residues for optimization during model building.
        
        Args:
            residues (tuple): Indices of mutated residues to optimize
            excluded (tuple): Indices of residues to keep fixed
            radius (float): Radius around mutations to include in optimization
        """
        self._opt_residues = tuple(map(int, residues))
        self._opt_residues_excluded = tuple(map(int, excluded))
        self._opt_residues_radius = radius

    def select_atoms(self):
        """
        Define atom selection for optimization.
        
        Selects atoms for refinement based on mutation positions,
        optimization radius, and exclusion rules. This overrides
        MODELLER's default behavior to optimize the entire structure.
        
        Returns:
            modeller.selection: Atom selection for optimization
        """
        # Start with mutated residues
        if self._opt_residues:
            residues = [self.residues[pos] for pos in self._opt_residues]
        else:
            residues = self.residues

        selection = mdl.selection(residues)

        # Expand selection by radius if specified
        if self._opt_residues_radius > 0:
            selection = selection.select_sphere(self._opt_residues_radius).by_residue()

        # Remove explicitly excluded residues
        excluded_residues = [self.residues[pos] for pos in self._opt_residues_excluded]
        excluded_selection = mdl.selection(excluded_residues).by_residue()

        return selection - excluded_selection
