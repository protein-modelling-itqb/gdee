"""
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
    def __init__(self, parameters):
        self.parameters = parameters
        mdl.log.none()
        self.env = mdl.environ()
        self.env.io.hetatm = True
        self.env.io.water = True
        self.env.edat.dynamic_lennard = True
        self.env.schedule_scale[mdl.physical.lennard_jones] = 1.0
        mdl.log.level(0, 0, 0, 0, 0)

    def run(self, job_data):
        temp_dir = TemporaryDirectory(prefix="gdee_modeller")
        temp_path = Path(temp_dir.name)
        Path(self.parameters["pdb_file"]).copy(temp_path / "template.pdb")

        with temp_path:
            self.write_alignment(job_data)

            model_data = self.build_models(job_data)
            if not model_data:
                job_data.fatal_error = True
                return job_data

            model_data.sort(key=lambda x: x["Normalized DOPE score"])
            top_models = model_data[:self.parameters["num_models"]]

            pdb_list = [str(model["name"]) for model in top_models]
            structure = mda.Universe(pdb_list[0], pdb_list)
            self.rename_models(structure, job_data.variant)

            model_list = []
            for i, ts in enumerate(structure.trajectory):
                filename = "model_{:04d}.pdb".format(i)
                structure.atoms.write(str(job_data.job_dir / filename))

                model = DataContainer()
                model.scores = DataContainer()
                model.scores.norm_dope = top_models[i]["Normalized DOPE score"]
                model.scores.molpdf = top_models[i]["molpdf"]
                model.pdb = filename
                model.rejected = False
                model_list.append(model)

        if model_list:
            modeled = DataContainer()
            modeled.method = "modeller"
            modeled.models = model_list
            job_data.modeling = modeled

        else:
            job_data.fatal_error = True

        return job_data

    def write_alignment(self, job_data):
        template = ">P1;template\nstructureX:template.pdb:.:.:.:.::::\n{}*\n\n>P1;model\nsequence:model.pdb:.:.:.:.::::\n{}*\n"
        with open("alignment.ali", "w") as fd:
            fd.write(template.format(
                job_data.wildtype.to_modeller(),
                job_data.variant.to_modeller()
            ))

    def build_models(self, job_data):
        model = MutationModel(
            self.env,
            alnfile="alignment.ali",
            knowns="template",
            sequence="model",
            assess_methods=(automodel.assess.normalized_dope),
        )

        model.starting_model = 1
        model.ending_model = 3 * self.parameters["num_models"]
        model.final_malign3d = True

        opt_level = self.parameters["optimize_level"]
        if opt_level == 0:
            model.library_schedule = automodel.autosched.very_fast
            model.md_level = automodel.refine.fast
            model.max_var_iterations = 100
            model.repeat_optimization = 1

        elif opt_level == 1:
            model.library_schedule = automodel.autosched.normal
            model.md_level = automodel.refine.slow
            model.max_var_iterations = 200
            model.repeat_optimization = 1

        else:
            model.library_schedule = automodel.autosched.slow
            model.md_level = automodel.refine.very_slow
            model.max_var_iterations = 300
            model.repeat_optimization = 1


        with contextlib.redirect_stdout(None):
            model.select_opt_residues(
                job_data.mut_index,
                job_data.fixed_index,
                self.parameters["optimize_radius"]
            )

            model.make()

        model_data = []
        for data in model.outputs:
            if "Normalized DOPE score" not in data:
                continue

            data["name"] = Path(data["name"]).stem + "_fit.pdb"
            model_data.append(data)

        return model_data

    def rename_models(self, structure, prot_seq):
        res_idx = 0
        for chain in prot_seq:
            segment = structure.add_Segment(segid=chain.code)

            for seq_pos in chain:
                if seq_pos.is_gap:
                    continue

                residue = structure.residues[res_idx]
                residue.resid = seq_pos.resid
                residue.resname = seq_pos.resname
                residue.segment = segment
                res_idx += 1


class MutationModel(automodel.automodel):
    def select_opt_residues(self, residues, excluded, coff):
        self._opt_residues = tuple(map(int, residues))
        self._opt_residues_excluded = tuple(map(int, excluded))
        self._opt_residues_coff = coff

    def select_atoms(self):
        if self._opt_residues:
            res = [self.residues[pos] for pos in self._opt_residues]
        else:
            res = self.residues

        sel = mdl.selection(res)

        if self._opt_residues_coff > 0:
            sel = sel.select_sphere(self._opt_residues_coff).by_residue()

        excluded_res = [self.residues[pos] for pos in self._opt_residues_excluded]
        excluded = mdl.selection(excluded_res).by_residue()

        return sel - excluded
