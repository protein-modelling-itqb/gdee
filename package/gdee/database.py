"""SQLite database interface for GDEE results storage."""

import sqlite3 as sql
import os
import uuid


def list_serialize(values):
    """Serialize list of values to pipe-delimited string.

    Args:
        values: List of numeric or string values

    Returns:
        str: Pipe-delimited serialized string
    """
    try:
        return "|".join(map(lambda x: "{:.4f}".format(x), values))
    except ValueError:
        pass

    return "|".join(map(lambda x: "{}".format(x), values))


class Database:
    """SQLite database interface for storing GDEE results.

    Manages storage of variants, models, docking results, and measurements
    with proper relationships and constraints.
    """

    def __init__(self, filename):
        if os.path.splitext(filename)[1] != ".sqlite3":
            filename += ".sqlite3"

        self.filename = filename
        self._conn = None
        self._metric_ids = {}

    def connect(self):
        # Connect to database
        self._conn = sql.connect(self.filename, timeout=120)
        self._conn.execute("PRAGMA foreign_keys = ON")

    def create_tables(self):
        """Create database schema with all required tables.

        Creates tables for proteins, variants, models, evaluations,
        poses, metrics, and measurements with proper relationships.
        """
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
        """Register a protein target.

        Args:
            name: Protein name
            uniprot: UniProt accession code (optional)

        Returns:
            int: Protein database ID
        """
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
        """Fetch all variants for a protein.

        Args:
            prot_id: Protein database ID

        Returns:
            list: List of (name, variant_id) tuples
        """
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
        """Check if variant exists.

        Args:
            prot_id: Protein database ID
            mutations: Variant name

        Returns:
            bool: True if variant exists
        """
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
        """Register a new variant.

        Args:
            prot_id: Parent protein database ID
            name: Variant name
            sequence: MODELLER format sequence
            directory: Working directory path
            wildtype: Whether variant is wildtype
            pdb_file: Template PDB file (optional)
            pdb_code: PDB code (optional)

        Returns:
            int: Variant database ID
        """
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
        """Register a 3D model.

        Args:
            variant_id: Parent variant database ID
            method: Modeling method name
            scores: JSON-formatted quality scores
            pdb_file: Path to model PDB file
            rejected: Whether model failed quality assessment

        Returns:
            int: Model database ID
        """
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
        """Register a docking evaluation.

        Args:
            variant_id: Variant database ID
            model_id: Model database ID
            evaluation: Evaluation data container

        Returns:
            int: Evaluation database ID
        """
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
        """Register docking poses with energies.

        Args:
            eval_id: Evaluation database ID
            energy: List of binding energies

        Returns:
            list: List of pose database IDs
        """
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
        """Register a measurement metric.

        Args:
            name: Metric name
            identifier: Unique metric identifier

        Returns:
            int: Metric database ID
        """
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

    def register_measurements(self, eval_id, pose_id_list, measurements):
        """Register computed measurements.

        Args:
            eval_id: Evaluation database ID
            pose_id_list: List of pose database IDs
            measurements: List of measurement data containers
        """
        cursor = self.conn.cursor()
        for measurer in measurements:
            metric_id = self.register_metric(measurer.name, measurer.identifier)

            for pose_id, value in zip(pose_id_list, measurer.data):
                cursor.execute(
                    "INSERT INTO"
                    "    Measurements ("
                    "        metric_id,"
                    "        pose_id,"
                    "        value"
                    "    ) "
                    "VALUES (?, ?, ?);",
                    (metric_id, pose_id, float(value))
                )

        self.conn.commit()


class Ranked_Database:
    def __init__(self, filename):
        if os.path.splitext(filename)[1] != ".sqlite3":
            filename += ".sqlite3"

        self.filename = filename
        self.table = None
        self._conn = None

    def connect(self):
        # Connect to database
        self._conn = sql.connect(self.filename, timeout=120)
        self._conn.execute("PRAGMA foreign_keys = ON")

    @property
    def conn(self):
        # Database is connected only when needed to allow
        # for multiple instantiations on MPI platform
        if self._conn is None:
            self.connect()
            if self.table is not None:
                self.create_table()

        return self._conn

    def fetch_ranked_variants(self, table):
        """
        Fetch all variants from ranked table.
        """
        cursor = self.conn.execute(
                "SELECT"
                "     * "
                "FROM"
                "     {}".format(table)
        )
        return cursor.fetchall()

    def _generate_table_name(self):
        """
        Generate a unique table name for temporary ranked results.
        """
        return "R" + str(uuid.uuid4()).replace("-", "")

    def create_table(self):
        # Create table
        with self.conn:
            self.conn.executescript(
                "CREATE TABLE IF NOT EXISTS"
                "    {} ("
                "        energy FLOAT NOT NULL,"
                "        name TEXT NOT NULL,"
                "        directory TEXT NOT NULL,"
                "        is_wildtype INT NOT NULL,"
                "        model_id INT NOT NULL,"
                "        variant_id INT NOT NULL,"
                "        eval_id INT NOT NULL,"
                "        pose_index INT NOT NULL,"
                "        docking_file TEXT NOT NULL"
                "    );".format(self.table)
        )

    def register_variant(self, energy, name, directory, wildtype, model_id, variant_id, eval_id, pose_index, docking_file):
        conn = self.conn
        cursor = conn.execute(
            "INSERT INTO"
            "    {} ("
            "        energy,"
            "        name,"
            "        directory,"
            "        is_wildtype,"
            "        model_id,"
            "        variant_id,"
            "        eval_id,"
            "        pose_index,"
            "        docking_file"
            "    ) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);".format(self.table),
            (energy, name, directory, bool(wildtype), model_id, variant_id, eval_id, pose_index, docking_file)
        )
        conn.commit()

    def by_num_mutations(self, num_mutations, table, table_name):
        like = "%" + "|%" * (num_mutations - 1)
        not_like = "%" + "|%" * num_mutations
        self._temp_table = self._generate_table_name()
        self.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT * "
            "FROM {} "
            "WHERE name LIKE '{}' AND name NOT LIKE '{}';".format(self._temp_table, table, like, not_like)
        )
        return self.export_sqlite(self.filename, table_name)

    def export_sqlite(self, file_name, table_name):
        self.conn.executescript(
                "ATTACH DATABASE '{0}' AS exportdb; "
                "DROP TABLE IF EXISTS exportdb.{1}; "
                "CREATE TABLE 'exportdb'.'{1}' AS "
                "SELECT * "
                "FROM {2} "
                "ORDER BY energy ASC; "
                "DETACH DATABASE exportdb;".format(file_name, table_name, self._temp_table)
        )
        return Ranked_Database(file_name)
