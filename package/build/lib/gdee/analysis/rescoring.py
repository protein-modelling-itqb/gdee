import uuid
import os
from oddt.virtualscreening import virtualscreening as vs
import math
from .filters import Rank
from ..database import Ranked_Database
from parallelbar import progress_map

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

def _generate_table_name():
    return "R" + str(uuid.uuid4()).replace("-", "")

class Rescore:
    def __init__(self, table, pickle_path, path, n_cpu=-1):
        self._table = table
        self._path = path
        self._pickle_path = pickle_path
        self._n_cpu = n_cpu

    def run(self, database):
        self.variants = list(database.fetch_ranked_variants(self._table))
        results = self.run_mp(self.variants)

        self._database = database
        self._temp_table = self.create_table()

        for result in results:
            self.update_table(result)

    def run_mp(self, arg, dlb=False):
        if self._n_cpu != -1:
            cores_available = self._n_cpu
        else:
            cores_available = None
            
        chunksize=max(len(arg) // cores_available, 1)
        if dlb: # Remove ?
            chunksize = max(chunksize/10, 1)
        
        scores = progress_map(self.run_variant, arg, n_cpu=cores_available, chunk_size=chunksize)

        return scores


    def run_variant(self, var):
        variant = Variant(var[0], var[1], var[2], var[3], var[4], self._path)
        score = variant.run_metamodel(self._pickle_path)
        
        return (variant.name, score)


    def create_table(self):
        table = _generate_table_name()
        self._database.conn.execute(
                "CREATE TEMP TABLE {} AS "
                "SELECT * "
                "FROM {}; ".format(table, self._table)
        )

        return table


    def update_table(self, values):
        self._database.conn.execute(
                "UPDATE {} "
                "SET energy = ? "
                "WHERE name = ?;".format(self._temp_table),
                (values[1], values[0])
        )


    def export_sqlite(self, file_name, table_name):
        self._database.conn.executescript(
                "ATTACH DATABASE '{0}' AS exportdb; "
                "DROP TABLE IF EXISTS exportdb.{1}; "
                "CREATE TABLE 'exportdb'.'{1}' AS "
                "SELECT * "
                "FROM {2} "
                "ORDER BY energy ASC; "
                "DETACH DATABASE exportdb;".format(file_name, table_name, self._temp_table)
        )


        return Rank(table_name, Ranked_Database(file_name))

            


class Variant:
    def __init__(self, energy, name, directory, pose, file, path):
        self._vina_score = energy
        self._name = name
        self._pose = pose
        self._docking_file = os.path.join(path, directory, file)
        self._model_file = os.path.join(path, directory, "model_{}".format(file[-8:]))

    @property
    def name(self):
        return self._name

    def kd_to_energy(self, pkd):
        R = 8.314
        T = 298.15
        kd = pow(10, -(float(pkd)))
        j = R * T * (math.log(float(kd)))
        kj = j / 1000
        kcal = kj / 4.18
        return kcal


    # TODO make training automated -> Goal: The user dont have to choose models 
    # TODO try to rescore only the correct pose
    def run_ml_functions(self, functions):
        ml_scores = []
        for function in functions:
            vs_rescore = vs(n_cpu=1)
            vs_rescore.load_ligands("pdb", self._docking_file)
            vs_rescore.score(function=function, protein=self._model_file)

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

            ml_scores.append(self.kd_to_energy(results[self._pose]))
        
        return ml_scores

    
    # TODO implement Vinardo base model ?
    def run_metamodel(self, base_functions):
        ml_scores = self.run_ml_functions(base_functions)

        w0 = 3.439315083935229
        wrf1 = -0.97639133
        wrf2 = 2.64010434
        wrf3 = -0.28160126
        wpnn = 0.06016545
        wvina = -0.04383675

        metamodel_score = w0 + (wrf1 * ml_scores[0]) + (wrf2 * ml_scores[1]) + (wrf3 * ml_scores[2]) + (wpnn * ml_scores[3]) + (wvina * self._vina_score)

        return metamodel_score

