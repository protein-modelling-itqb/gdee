import os
from oddt.virtualscreening import virtualscreening as vs
import math
from ..misc import DataContainer
from ..database import Ranked_Database
import sys
import traceback

import warnings
warnings.filterwarnings("ignore", category=UserWarning)


class Rescore:
    def __init__(self, output_db, table):
        self.table = table
        self.path = None
        self.pickle_path = None
        self.output_db = Ranked_Database(output_db)
        self.output_db.table = self.table

    def get_variants(self, db):
        variant_list = list(db.fetch_ranked_variants(self.table))
        self.variant_it = iter(variant_list)

    def fetch_next_job(self):
        try:
            variant = next(self.variant_it)
        except StopIteration:
            return None

        energy = variant[0]
        name = variant[1]
        variant_dir = self.path / variant[2]
        wt = variant[3]
        model_id = variant[4]
        variant_id = variant[5]
        eval_id = variant[6]
        pose = variant[7]
        docking_file = variant[8]

        job = DataContainer()
        job.variant_dir = variant_dir
        job.energy = energy
        job.name = name
        job.pose = pose
        job.docking_file = docking_file
        job.wt = wt
        job.model_id = model_id
        job.variant_id = variant_id
        job.eval_id = eval_id
        return job

    def run(self, job_data):
        docking_file = self.path / job_data.variant_dir / job_data.docking_file
        model_file = self.path / job_data.variant_dir / "model_{}".format(job_data.docking_file[-8:])
        pose = job_data.pose
        vina_score = job_data.energy
        name = job_data.name

        results = self.run_metamodel(docking_file, model_file, pose, vina_score, name)

        job_data.results = results
        return job_data

    def kd_to_energy(self, pkd):
        R = 8.314
        T = 298.15
        kd = pow(10, -(float(pkd)))
        j = R * T * (math.log(float(kd)))
        kj = j / 1000
        kcal = kj / 4.18
        return kcal

    def run_ml_functions(self, docking_file, model_file, pose, name):
        ml_scores = []
        for function in self.pickle_path:
            vs_rescore = vs(n_cpu=1)
            vs_rescore.load_ligands("pdb", docking_file)
            vs_rescore.score(function=function, protein=model_file)

            # Save result
            results = []
            for mol in vs_rescore.fetch():
                data = mol.data.to_dict()

                if len(data) > 0:
                    data['name'] = mol.title
                else:
                    print('There is no data', file=sys.stderr)
                    return False

                for key in data:
                    if key.startswith('rf'):
                        results.append(data[key])
                    elif key.startswith('PLEC'):
                        results.append(data[key])

            try:
                ml_scores.append(self.kd_to_energy(results[pose]))

            except Exception as error:
                print("Exception caught for variant {}:".format(name), traceback.print_exc(), error, sep="\n")
        return ml_scores

    def run_metamodel(self, docking_file, model_file, pose, vina_score, name):
        ml_scores = self.run_ml_functions(docking_file, model_file, pose, name)

        w0 = 3.439315083935229
        wrf1 = -0.97639133
        wrf2 = 2.64010434
        wrf3 = -0.28160126
        wpnn = 0.06016545
        wvina = -0.04383675

        metamodel_score = w0 + (wrf1 * ml_scores[0]) + (wrf2 * ml_scores[1]) + (wrf3 * ml_scores[2]) + (wpnn * ml_scores[3]) + (wvina * vina_score)
        return metamodel_score

    def save_results(self, data):
        name = data.name
        variant_dir = os.path.basename(data.variant_dir)

        self.output_db.register_variant(data.results, name, variant_dir, data.wt, data.model_id, data.variant_id, data.eval_id, data.pose, data.docking_file)

        print("Ended with variant '{}'".format(name))
        return True
