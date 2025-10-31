from gdee.measurement import MeasurerFactory
from gdee.engineer import Ligand
from gdee.database import Database
import multiprocessing as mp
from path import Path


class PostPipelineMeasurer:
    def __init__(self):
        self.metrics = None

    def measure_saved_results(self, variant_dir, name):
        """
        Perform measurements on saved results after the pipeline has finished.

        Args:
            variant_dir (Path): Directory containing the variant's results.
            name (str): Name of the ligand.

        Returns:
            dict: A dictionary of measurement results, keyed by (variant_dir, protein_file, ligand_file).
        """
        # Find all protein and ligand files
        protein_files = list(variant_dir.glob("model_*.pdb"))
        ligand_files = list(variant_dir.glob("docking_*.pdb"))
        if not protein_files or not ligand_files:
            print(f"Skipping {variant_dir}, missing required files.")
            return {}

        # Map four-digit keys to protein and ligand files
        protein_map = {f.stem.split("_")[-1]: f for f in protein_files}
        ligand_map = {}
        for f in ligand_files:
            key = f.stem.split("_")[-1]
            ligand_map.setdefault(key, []).append(f)

        # Find matching keys
        matching_keys = protein_map.keys() & ligand_map.keys()
        if not matching_keys:
            print(f"Skipping {variant_dir}, no matching protein and ligand files.")
            return {}

        # Perform measurements for each matching pair
        results = {}
        for key in matching_keys:
            protein_file = protein_map[key]
            for ligand_file in ligand_map[key]:
                try:
                    # Create a ligand object and add measurements
                    ligand = Ligand(name=name, filename=ligand_file)
                    for metric in self.metrics[name]:
                        ligand.add_measurement(metric[0], metric[1], metric[2], metric[3])

                    # Create a measurer and perform measurements
                    measurer_factory = MeasurerFactory()
                    measurer_factory.ligand = ligand
                    measurer = measurer_factory.make()
                    measurements = measurer.measure_from_files(protein_file, ligand_file)

                    if not measurements:
                        print(f"No measurements found for {variant_dir} (key {key}, ligand {ligand_file.name})")
                        continue

                    results[(variant_dir.name, protein_file.name, ligand_file.name)] = measurements
                except Exception as error:
                    print(f"Error measuring {variant_dir} (key {key}, ligand {ligand_file.name}): {error}")
        
        return results

    def run_measurements(self, results_dir, name, database_path):
        """
        Run measurements on saved results and register them in the database.

        Args:
            results_dir (str): Path to the directory containing saved results.
            name (str): Name of the ligand.
            database_path (str): Path to the database file.
        """
        database = Database(database_path)
        results_dir = Path(results_dir)

        # Collect all variant directories
        jobs = [variant_dir for variant_dir in results_dir.iterdir() if variant_dir.is_dir()]

        print("Starting measurements for all variants...")
        # Use multiprocessing to measure results in parallel
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = pool.starmap(self.measure_saved_results, [(variant_dir, name) for variant_dir in jobs])
                    
        print("Finished measuring all variants.")
        # Register the measurements in the database
        for result in results:
            if result:
                for key, measurements in result.items():
                    try:
                        database.register_new_measurements(measurements, key[0], key[2])
                    except Exception as error:
                        print(f"Error registering measurements for {key}: {error}")
    
        # Stop the worker thread after all results are processed
        database.data_queue.join()
        database.stop_worker()

