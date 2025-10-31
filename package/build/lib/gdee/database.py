"""
"""


import sqlite3 as sql
import os
from queue import Queue
import threading


def list_serialize(values):
    try:
        return "|".join(map(lambda x: "{:.4f}".format(x), values))
    except ValueError:
        pass

    return "|".join(map(lambda x: "{}".format(x), values))


class Database:
    def __init__(self, filename):
        if os.path.splitext(filename)[1] != ".sqlite3":
            filename += ".sqlite3"

        self.filename = filename
        self._conn = None
        self._metric_ids = {}
        self.data_queue = Queue()
        self.worker_thread = threading.Thread(target=self._worker_thread, daemon=True)
        self.worker_thread.start()

    def connect(self):
        # Connect to database
        self._conn = sql.connect(self.filename, timeout=120)
        self._conn.execute("PRAGMA foreign_keys = ON")

    def create_tables(self):
        # Create tables
        with self.conn:
            self.conn.executescript(
                "CREATE TABLE IF NOT EXISTS"
                "    Proteins ("
                "        prot_id INTEGER PRIMARY KEY,"
                "        uniprot TEXT,"
                "        name TEXT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    ProteinMetadata ("
                "        protmeta_id INTEGER PRIMARY KEY,"
                "        prot_id"
                "            REFERENCES Proteins(prot_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        data TEXT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Variants ("
                "        variant_id INTEGER PRIMARY KEY,"
                "        name TEXT NOT NULL,"
                "        sequence TEXT NOT NULL,"
                "        prot_id"
                "            REFERENCES Proteins(prot_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        directory TEXT,"
                "        is_wildtype INT NOT NULL,"
                "        pdb_file TEXT,"
                "        pdb_code TEXT,"
                "        UNIQUE("
                "            name,"
                "            sequence,"
                "            prot_id"
                "        )"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Models ("
                "        model_id INTEGER PRIMARY KEY,"
                "        variant_id"
                "            REFERENCES Variants(variant_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        method TEXT NOT NULL,"
                "        scores TEXT NOT NULL,"
                "        pdb_file TEXT NOT NULL,"
                "        rejected INT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Evaluations ("
                "        eval_id INTEGER PRIMARY KEY,"
                "        variant_id"
                "            REFERENCES Variants(variant_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        model_id"
                "            REFERENCES Models(model_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        ligand_name TEXT NOT NULL,"
                "        ligand_file TEXT NOT NULL,"
                "        method TEXT NOT NULL,"
                "        pdb_file TEXT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Poses ("
                "        pose_id INTEGER PRIMARY KEY,"
                "        eval_id"
                "            REFERENCES Evaluations(eval_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        pdb_index INTEGER NOT NULL,"
                "        energy REAL NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Metrics ("
                "        metric_id INTEGER PRIMARY KEY,"
                "        name TEXT UNIQUE NOT NULL,"
                "        identifier TEXT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Measurements ("
                "        measurement_id INTEGER PRIMARY KEY,"
                "        metric_id"
                "            REFERENCES Metrics(metric_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        pose_id"
                "            REFERENCES Poses(pose_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        value REAL NOT NULL"
                "    );"
        )

    @property
    def conn(self):
        # Database is connected only when needed to allow
        # for multiple instantiations on MPI platform
        if self._conn is None:
            self.connect()
            self.create_tables()

        return self._conn

    def register_protein(self, name, uniprot=None):
        conn = self.conn
        cursor = conn.execute(
            "SELECT"
            "    prot_id "
            "FROM"
            "    Proteins "
            "WHERE"
            "    uniprot IS ?"
            "    AND"
            "    name IS ?;",
            (uniprot, name)
        )
        data = cursor.fetchone()

        if data:
            return data[0]

        cursor = conn.execute(
            "INSERT INTO"
            "    Proteins ("
            "        uniprot,"
            "        name"
            "    ) "
            "VALUES (?, ?);",
            (uniprot, name)
        )
        conn.commit()

        return cursor.lastrowid

    def fetch_variants(self, prot_id):
        cursor = self.conn.execute(
            "SELECT"
            "    name "
            "FROM"
            "    Variants "
            "WHERE"
            "    prot_id = ?;",
            (prot_id,)
        )
        return cursor.fetchall()

    def variant_exists(self, prot_id, mutations):
        cursor = self.conn.execute(
            "SELECT EXISTS ("
            "    SELECT"
            "        1"
            "    FROM"
            "        Variants"
            "    WHERE"
            "        prot_id = ?"
            "        AND"
            "        name = ?"
            ");",
            (prot_id, mutations,)
        )
        return bool(cursor.fetchone()[0])

    def register_variant(self, prot_id, name, sequence, directory, wildtype, pdb_file=None, pdb_code=None):
        conn = self.conn
        cursor = conn.execute(
            "INSERT INTO"
            "    Variants ("
            "        name,"
            "        sequence,"
            "        prot_id,"
            "        directory,"
            "        is_wildtype,"
            "        pdb_file,"
            "        pdb_code"
            "    ) "
            "VALUES (?, ?, ?, ?, ?, ?, ?);",
            (name, sequence, prot_id, directory, bool(wildtype), pdb_file, pdb_code)
        )
        conn.commit()

        return cursor.lastrowid

    def remove_variant(self, variant_name):
        conn = self.conn
        cursor = conn.execute(
            "DELETE FROM"
            "    Variants "
            "WHERE"
            "    name = ?;",
            (variant_name,)
        )
        conn.commit()

        return cursor.lastrowid

    def register_model(self, variant_id, method, scores, pdb_file, rejected):
        conn = self.conn
        cursor = conn.execute(
            "INSERT INTO"
            "    Models ("
            "        variant_id,"
            "        method,"
            "        scores,"
            "        pdb_file,"
            "        rejected"
            "    ) "
            "VALUES (?, ?, ?, ?, ?);",
            (variant_id, method, scores, pdb_file, rejected)
        )

        conn.commit()
        return cursor.lastrowid

    def register_evaluation(self, variant_id, model_id, evaluation):
        conn = self.conn
        cursor = conn.execute(
            "INSERT INTO"
            "    Evaluations ("
            "        variant_id,"
            "        model_id,"
            "        ligand_name,"
            "        ligand_file,"
            "        method,"
            "        pdb_file"
            "    ) "
            "VALUES (?, ?, ?, ?, ?, ?);",
            (variant_id, model_id, evaluation.ligand_name,
             evaluation.ligand_file, evaluation.method, evaluation.pdb)
        )
        return cursor.lastrowid

    def register_poses(self, eval_id, energy):
        pose_id = []
        conn = self.conn
        cursor = conn.cursor()
        for index, value in enumerate(energy):
            cursor.execute(
                "INSERT INTO"
                "    Poses ("
                "        eval_id,"
                "        pdb_index,"
                "        energy"
                "    ) "
                "VALUES (?, ?, ?);",
                (eval_id, index, value)
            )
            pose_id.append(cursor.lastrowid)

        conn.commit()
        return pose_id

    def fetch_metric_id(self, name):
        cursor = self.conn.execute(
            "SELECT"
            "    metric_id "
            "FROM"
            "    Metrics "
            "WHERE"
            "    name = ?;",
            (name,)
        )
        data = cursor.fetchall()
        if data:
            return data[0][0]

        return None

    def register_metric(self, name, identifier):
        if name not in self._metric_ids:
            conn = self.conn
            metric_id = self.fetch_metric_id(name)

            if metric_id is not None:
                self._metric_ids[name] = metric_id
            else:
                cursor = conn.execute(
                    "INSERT INTO"
                    "    Metrics ("
                    "        name,"
                    "        identifier"
                    "    ) "
                    "VALUES (?, ?);",
                    (name, identifier)
                )
                self._metric_ids[name] = cursor.lastrowid
                conn.commit()

        return self._metric_ids[name]


    def fetch_variant(self, directory):
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT variant_id FROM Variants WHERE directory = ?;",
            (directory,)
        )
        data = cursor.fetchone()
        return data[0] if data else None
    
    def fetch_evaluation(self, variant_id, ligand_file):
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT eval_id FROM Evaluations WHERE variant_id = ? AND pdb_file = ?;",
            (variant_id, ligand_file)
        )
        data = cursor.fetchone()
        return data[0] if data else None
    
    def fetch_pose(self, eval_id, pose_index):
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT pose_id FROM Poses WHERE eval_id = ? AND pdb_index = ?;",
            (eval_id, pose_index)
        )
        data = cursor.fetchone()
        return data[0] if data else None
        
    def register_new_measurements(self, measurements, variant_dir, ligand_file):
        """
        Add new measurements to the queue for asynchronous database writing.
        """
        self.data_queue.put((measurements, variant_dir, ligand_file))

    def _worker_thread(self):
        """
        Thread that processes write requests from the queue and writes them to the database.
        """
        self.connect()
        self._conn.execute("PRAGMA journal_mode = WAL") 
        cursor = self.conn.cursor()

        while True:
            data = self.data_queue.get()
            if data is None:  # Sentinel value to terminate
                break

            try:
                self._write_measurements(cursor, *data)
            except sql.Error as e:
                print(f"Error writing to database: {e}")

            self.data_queue.task_done()

        self.conn.close()

    def _write_measurements(self, cursor, measurements, variant_dir, ligand_file):
        """
        Register new measurements in the database in batches.
        """
        cursor = self.conn.cursor()
        variant_id = self.fetch_variant(variant_dir)
        if variant_id is None:
            print(f"Variant not found for {variant_dir}")
            return
        eval_id = self.fetch_evaluation(variant_id, ligand_file)
        if eval_id is None:
            print(f"Evaluation not found for {variant_dir} and {ligand_file}")
            return

        batch_data = []
        for measurer in measurements:
            metric_id = self.register_metric(measurer.name, measurer.identifier)
            for pdb_index, value in enumerate(measurer.data):
                pose_id = self.fetch_pose(eval_id, pdb_index)
                batch_data.append((metric_id, pose_id, float(value)))

        # Perform batch insert
        cursor.executemany(
            "INSERT INTO Measurements (metric_id, pose_id, value) VALUES (?, ?, ?);",
            batch_data
        )

        print(f"Inserted {len(batch_data)} measurements")
        self.conn.commit()

    def stop_worker(self):
        """
        Stop the worker thread and wait for it to finish.
        """
        self.data_queue.put(None)
        self.worker_thread.join()