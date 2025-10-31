"""
"""


import numpy as np


class PDBQT:
    def __init__(self, filename):
        self.atom_data = []
        self.atom_qt = []
        self.models = []
        self.parse(filename)

    def size(self):
        return len(self.models)

    def __getitem__(self, key):
        return self.models[key]

    def __iter__(self):
        return iter(self.models)

    def parse(self, filename):
        coords = []
        charges = []
        types = []
        energy = 0
        has_models = False
        first_pass = True

        with open(filename) as file:
            for line in file:
                trimmed = line.strip()

                if trimmed.startswith("MODEL"):
                    has_models = True
                    coords.clear()
                    charges.clear()
                    types.clear()

                elif trimmed.startswith("USER") and \
                    "Estimated Free Energy of Binding" in trimmed:
                    columns = trimmed.split()
                    energy = float(columns[7])

                elif trimmed.startswith("REMARK"):
                    columns = trimmed.split()
                    if "VINA" in columns and "RESULT:" in columns:
                        energy = float(columns[3])

                    elif "minimizedAffinity" in columns:
                        energy = float(columns[2])

                elif trimmed.startswith("ATOM") or trimmed.startswith("HETATM"):
                    if first_pass:
                        self.atom_data.append(trimmed[:27])
                        self.atom_qt.append(trimmed[70:])

                    x = float(trimmed[30:38])
                    y = float(trimmed[38:46])
                    z = float(trimmed[46:54])
                    coords.append((x, y, z))

                elif trimmed.startswith("ENDMDL"):
                    first_pass = False
                    if coords:
                        self.models.append(DockingModel(coords, energy))

            if not has_models:
                if coords:
                    self.models.append(DockingModel(coords, energy))

    def write_pdb(self, filename):
        with open(filename, "w") as fd:
            for model_idx, model in enumerate(self.models):
                fd.write("MODEL {:d}\nREMARK ENERGY: {:.2f}\n".format(model_idx, model.energy))

                for atom_idx in range(len(self.atom_data)):
                    fd.write("{:<30s}{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
                        self.atom_data[atom_idx],
                        model.coords[atom_idx, 0],
                        model.coords[atom_idx, 1],
                        model.coords[atom_idx, 2],
                        0,
                        model.energy,
                    ))

                fd.write("ENDMDL\n")


class DockingModel:
    def __init__(self, coords, energy):
        self.coords = np.array(coords, "float64")
        self.energy = energy

    def rmsd(self, other):
        if self.coords.shape != other.coords.shape:
            raise ValueError("Models must have the same number of atoms")

        return np.sqrt(np.sum((self.coords - other.coords) ** 2) / self.coords.shape[0])

    def centroid(self):
        return self.coords.mean(0)

    def box(self):
        return self.coords.max(0) - self.coords.min(0)
