"""
SQLite database interface for GDEE platform.

This module provides database operations for storing and retrieving protein variants,
3D models, docking evaluations, poses, and measurement data from the GDEE platform.
"""


import sqlite3 as sql
import os


def list_serialize(values):
    """
    Serialize a list of values to a pipe-delimited string format.
    
    Attempts to format numeric values to 4 decimal places, falls back
    to string representation for non-numeric values.
    
    Args:
        values (list): List of values to serialize
        
    Returns:
        str: Pipe-delimited string representation of values
    """
    try:
        return "|".join(map(lambda x: "{:.4f}".format(x), values))
    except ValueError:
        pass

    return "|".join(map(lambda x: "{}".format(x), values))


class Database:
    """
    SQLite database interface for GDEE platform data.

    This class manages the storage and retrieval of GDEE platform
    results including variants, 3D models, docking evaluations, poses, and measurements.
    The database schema supports hierarchical data relationships and foreign key constraints.
    
    Database Schema:
        - Proteins: Target protein information
        - ProteinMetadata: Additional protein metadata
        - Variants: Protein sequence variants
        - Models: 3D structural models (MODELLER output)
        - Evaluations: Molecular docking results
        - Poses: Individual docking poses with energies
        - Metrics: Metric definitions
        - Measurements: Distance measurements
    
    Attributes:
        filename (str): Path to SQLite database file
        _conn: SQLite connection object (lazy initialized)
        _metric_ids (dict): Cache for metric ID lookups
    """
    
    def __init__(self, filename):
        """
        Initialize database connection manager.
        
        Args:
            filename (str): Path to SQLite database file (adds .sqlite3 extension if missing)
        """
        if os.path.splitext(filename)[1] != ".sqlite3":
            filename += ".sqlite3"

        self.filename = filename
        self._conn = None
        self._metric_ids = {}

    def connect(self):
        """
        Establish connection to SQLite database.
        
        Enables foreign key constraints and sets a 120-second timeout for database locks.
        """
        self._conn = sql.connect(self.filename, timeout=120)
        self._conn.execute("PRAGMA foreign_keys = ON")

    def create_tables(self):
        """
        Create the complete database schema for GDEE platform.
        
        Creates tables with foreign key relationships:
        - Proteins → Variants → Models → Evaluations → Poses
        - Metrics → Measurements (linked to Poses)
        
        All tables use cascading deletes and updates for data integrity.
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
        """
        Get database connection with lazy initialization.
        
        Database connection and table creation are deferred until first access
        to support MPI platform with multiple process instantiation.
        
        Returns:
            sqlite3.Connection: Active database connection
        """
        if self._conn is None:
            self.connect()
            self.create_tables()

        return self._conn

    def register_protein(self, name, uniprot=None):
        """
        Register a target protein in the database.
        
        Args:
            name (str): Protein name identifier
            uniprot (str, optional): UniProt accession number
            
        Returns:
            int: Database protein ID (prot_id)
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
        """
        Retrieve all variant names for a given protein.
        
        Args:
            prot_id (int): Database protein ID
            
        Returns:
            list: List of tuples containing variant names
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
        """
        Check if a variant already exists for a protein.
        
        Args:
            prot_id (int): Database protein ID
            mutations (str): Variant name/mutation string
            
        Returns:
            bool: True if variant exists, False otherwise
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
        """
        Register a protein variant in the database.
        
        Args:
            prot_id (int): Database protein ID
            name (str): Variant name/mutation identifier
            sequence (str): Protein sequence
            directory (str): Working directory path
            wildtype (bool): Whether this is the wildtype sequence
            pdb_file (str, optional): PDB template file path
            pdb_code (str, optional): PDB accession code
            
        Returns:
            int: Database variant ID (variant_id)
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
        """
        Remove a variant and all associated data from the database.
        
        Args:
            variant_name (str): Name of variant to remove
            
        Returns:
            int: Last row ID from delete operation
        """
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
        """
        Register a 3D structural model (from ModellerBuilder).
        
        Args:
            variant_id (int): Database variant ID
            method (str): Modeling method (e.g., "modeller")
            scores (str): JSON string of quality scores (DOPE, VoroMQA)
            pdb_file (str): Path to model PDB file
            rejected (bool): Whether model failed quality assessment
            
        Returns:
            int: Database model ID (model_id)
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
        """
        Register a molecular docking evaluation.
        
        Args:
            variant_id (int): Database variant ID
            model_id (int): Database model ID
            evaluation: Evaluation data container with docking results
            
        Returns:
            int: Database evaluation ID (eval_id)
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
        """
        Register docking poses with their binding energies.
        
        Args:
            eval_id (int): Database evaluation ID
            energy (list): List of binding energies for each pose
            
        Returns:
            list: List of database pose IDs (pose_id)
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
        """
        Retrieve database ID for a measurement metric.
        
        Args:
            name (str): Metric name
            
        Returns:
            int or None: Database metric ID if found, None otherwise
        """
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
        """
        Register a measurement metric type (e.g., distance calculations).
        
        Uses caching to avoid repeated database queries for the same metric.
        
        Args:
            name (str): Metric name
            identifier (str): Metric identifier string (includes selection strings)
            
        Returns:
            int: Database metric ID (metric_id)
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
        """
        Register measurement values for docking poses.
        
        Stores distance calculations and other geometric measurements computed
        by the Measurer component for pose filtering and analysis.
        
        Args:
            eval_id (int): Database evaluation ID
            pose_id_list (list): List of database pose IDs
            measurements (list): List of measurement data containers
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
