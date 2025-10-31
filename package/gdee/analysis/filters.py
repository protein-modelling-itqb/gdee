"""
Advanced filtering and ranking system for GDEE result analysis.

This module provides a powerful SQL-based filtering system for analyzing protein
engineering results. It supports complex queries combining multiple measurements,
variant properties, and docking energies with boolean logic operations.
"""

import csv
import uuid
import numbers
from ..database import Database


def _generate_table_name():
    """
    Generate unique temporary table name.
    
    Creates a unique identifier for temporary SQL tables used in
    filtering operations to avoid naming conflicts.
    
    Returns:
        str: Unique table name starting with 'R'
    """
    return "R" + str(uuid.uuid4()).replace("-", "")


class Metric:
    """
    Measurement metric for filtering and comparison operations.
    
    This class represents a measurement type (e.g., distance calculations)
    and provides comparison operators for creating filtering rules. It enables
    complex queries like "distance < 5.0" or "metric1 > metric2".
    
    Attributes:
        _name (str): Metric name as stored in database
        _database (Database): Database connection for queries
        _metric_id (int): Database metric ID for efficient queries
    """
    
    def __init__(self, name, database):
        """
        Initialize metric for filtering operations.
        
        Args:
            name (str): Metric name (must exist in database)
            database (Database): Database connection
            
        Raises:
            RuntimeError: If metric is not found in database
        """
        self._name = name
        self._database = database
        self._metric_id = database.fetch_metric_id(name)
        if self._metric_id is None:
            raise RuntimeError("Metric '{}' not found".format(name))

    @property
    def name(self):
        """Get metric name."""
        return self._name

    def __lt__(self, value):
        """Less than comparison operator."""
        return self._compare("<", value)

    def __le__(self, value):
        """Less than or equal comparison operator."""
        return self._compare("<=", value)

    def __gt__(self, value):
        """Greater than comparison operator."""
        return self._compare(">", value)

    def __ge__(self, value):
        """Greater than or equal comparison operator."""
        return self._compare(">=", value)

    def __eq__(self, value):
        """Equality comparison operator."""
        return self._compare("=", value)

    def __ne__(self, value):
        """Not equal comparison operator."""
        return self._compare("!=", value)

    def _compare(self, operator, value):
        """
        Create comparison rule for metric filtering.
        
        Args:
            operator (str): SQL comparison operator
            value (number or Metric): Value or metric to compare against
            
        Returns:
            Rule: Filtering rule for pose selection
        """
        if isinstance(value, numbers.Number):
            return self._compare_number(operator, value)

        elif isinstance(value, Metric):
            return self._compare_metrics(operator, value)

    def _compare_number(self, operator, value):
        """
        Compare metric values against numeric threshold.
        
        Creates SQL query to find poses where metric satisfies
        the comparison condition with a numeric value.
        
        Args:
            operator (str): SQL comparison operator (<, >, =, etc.)
            value (number): Threshold value for comparison
            
        Returns:
            Rule: Filtering rule with matching pose IDs
        """
        table = _generate_table_name()
        self._database.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT pose_id "
            "FROM Measurements "
            "WHERE metric_id = ? AND value {} ?;".format(table, operator),
            (self._metric_id, value)
        )
        return Rule(table, self._database)

    def _compare_metrics(self, operator, other):
        """
        Compare two metrics against each other.
        
        Creates SQL query to find poses where one metric satisfies
        the comparison condition relative to another metric.
        
        Args:
            operator (str): SQL comparison operator
            other (Metric): Other metric for comparison
            
        Returns:
            Rule: Filtering rule with matching pose IDs
        """
        table = _generate_table_name()
        self._database.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT M1.pose_id "
            "FROM Measurements M1 "
            "INNER JOIN Measurements M2 "
            "ON M1.pose_id = M2.pose_id "
            "WHERE M1.metric_id = ? and M2.metric_id = ? "
            "AND M1.value {} M2.value;".format(table, operator),
            (self._metric_id, other._metric_id)
        )
        return Rule(table, self._database)


class Rule:
    """
    Filtering rule for pose selection and boolean operations.
    
    This class represents a set of pose IDs that satisfy filtering criteria.
    It supports boolean operations (AND, OR, NOT) to combine multiple
    filtering conditions into complex queries.
    
    Attributes:
        _table (str): Temporary table name containing pose IDs
        _database (Database): Database connection for operations
    """
    
    def __init__(self, table, database):
        """
        Initialize filtering rule.
        
        Args:
            table (str): Temporary table name with pose IDs
            database (Database): Database connection
        """
        self._table = table
        self._database = database

    def __and__(self, other):
        """
        Logical AND operation between rules.
        
        Creates new rule containing poses that satisfy both conditions.
        Equivalent to set intersection operation.
        
        Args:
            other (Rule): Other filtering rule
            
        Returns:
            Rule: Combined rule with intersection of pose IDs
        """
        table = _generate_table_name()
        self._database.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT pose_id "
            "FROM {} "
            "WHERE pose_id IN (SELECT pose_id FROM {});".format(table, self._table, other._table))
        return Rule(table, self._database)

    def __or__(self, other):
        """
        Logical OR operation between rules.
        
        Creates new rule containing poses that satisfy either condition.
        Equivalent to set union operation.
        
        Args:
            other (Rule): Other filtering rule
            
        Returns:
            Rule: Combined rule with union of pose IDs
        """
        table = _generate_table_name()
        self._database.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT pose_id "
            "FROM {} "
            "UNION "
            "SELECT pose_id "
            "FROM {} ".format(table, self._table, other._table))
        return Rule(table, self._database)

    def __invert__(self):
        """
        Logical NOT operation (complement).
        
        Creates new rule containing poses that do NOT satisfy the condition.
        Equivalent to set complement operation.
        
        Returns:
            Rule: Inverted rule with complementary pose IDs
        """
        table = _generate_table_name()
        self._database.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT pose_id "
            "FROM Measurements "
            "WHERE pose_id NOT IN (SELECT pose_id FROM {});".format(table, self._table))
        return Rule(table, self._database)

    def __bool__(self):
        """
        Prevent accidental boolean conversion.
        
        Rules should not be used in boolean contexts to avoid confusion
        with logical operations. Use bitwise operators instead.
        
        Raises:
            TypeError: Always raises to prevent misuse
        """
        raise TypeError("Rule not convertible to bool. Use the bitwise operators ~ (not), & (and), | (or) for boolean expressions.")

    def rank(self, ascending=True):
        """
        Convert rule to ranking with energy-based ordering.
        
        Creates a ranking of variants based on best docking energies
        for poses that satisfy the filtering rule.
        
        Args:
            ascending (bool): Sort order (True=lowest energy first)
            
        Returns:
            Rank: Ranked results ready for export or further analysis
        """
        return Rank.from_rule(self._table, self._database, ascending)


class Rank:
    """
    Ranked results for variant analysis and export.
    
    This class manages ranked lists of protein variants based on filtering
    criteria and energy scores. It provides methods for further filtering
    by mutation count, wildtype status, and data export capabilities.
    
    Attributes:
        _table (str): Temporary table name containing ranked results
        _database (Database): Database connection for operations
    """
    
    def __init__(self, table, database):
        """
        Initialize ranking system.
        
        Args:
            table (str): Temporary table name with ranked data
            database (Database): Database connection
        """
        self._table = table
        self._database = database

    @staticmethod
    def from_rule(poses_table, database, ascending):
        """
        Create ranking from filtering rule.
        
        This method:
        1. Joins poses with evaluation and variant data
        2. Groups by variant to find best energy per variant
        3. Orders results by energy (ascending/descending)
        4. Includes variant metadata for analysis
        
        Args:
            poses_table (str): Table name with filtered pose IDs
            database (Database): Database connection
            ascending (bool): Sort order for energy ranking
            
        Returns:
            Rank: Configured ranking ready for analysis
        """
        order = "ASC" if ascending else "DESC"
        temp_table = _generate_table_name()
        table = _generate_table_name()
        database.conn.executescript(
            "CREATE TEMP TABLE {0} AS "
            "    SELECT * "
            "    FROM Poses "
            "    WHERE pose_id IN (SELECT pose_id FROM {1}); "
            ""
            "CREATE TEMP TABLE {2} AS "
            "SELECT MIN(energy) energy, name, directory, is_wildtype, "
            "       model_id, Evaluations.variant_id, Evaluations.eval_id, "
            "       S.pdb_index pose_index, Evaluations.pdb_file docking_file "
            "FROM {0} S "
            "INNER JOIN Evaluations ON S.eval_id = Evaluations.eval_id "
            "INNER JOIN Variants ON Evaluations.variant_id = Variants.variant_id "
            "GROUP BY Evaluations.variant_id "
            "ORDER BY energy {3}; "
            ""
            "DROP TABLE IF EXISTS {0};".format(temp_table, poses_table, table, order)
        )
        return Rank(table, database)

    def by_num_mutations(self, num_mutations):
        """
        Filter ranking by number of mutations.
        
        Uses SQL LIKE patterns to match variant names with specific
        mutation counts based on the "|" separator in variant names.
        
        Args:
            num_mutations (int): Exact number of mutations to filter
            
        Returns:
            Rank: Filtered ranking with specified mutation count
            
        Examples:
            >>> rank.by_num_mutations(2)  # Only double mutants
            >>> rank.by_num_mutations(1)  # Only single mutants
        """
        like = "%" + "|%" * (num_mutations - 1)
        not_like = "%" + "|%" * num_mutations
        table = _generate_table_name()
        self._database.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT * "
            "FROM {} "
            "WHERE name LIKE '{}' AND name NOT LIKE '{}';".format(table, self._table, like, not_like)
        )
        return Rank(table, self._database)

    def by_wildtype(self):
        """
        Filter ranking to include only wildtype sequences.
        
        Selects variants marked as wildtype in the database.
        Useful for baseline comparisons and reference values.
        
        Returns:
            Rank: Filtered ranking containing only wildtype entries
        """
        table = _generate_table_name()
        self._database.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT * "
            "FROM {} "
            "WHERE is_wildtype = 1;".format(table, self._table)
        )
        return Rank(table, self._database)

    def export_csv(self, file_name, max_lines=100):
        """
        Export ranking to CSV file.
        
        Creates a CSV file with ranked results including energy scores,
        variant names, file paths, and metadata. Includes rank numbers
        for easy identification of top performers.
        
        Args:
            file_name (str): Output CSV file path
            max_lines (int or None): Maximum number of lines to export
                                    (None exports all results)
        """
        with open(file_name, "w") as fd:
            writer = csv.writer(fd)
            writer.writerow(("rank", "energy", "name", "directory", "is_wildtype", 
                           "model_id", "variant_id", "eval_id", "pose_index", "docking_file"))
            cursor = self._database.conn.execute(
                "SELECT ROW_NUMBER() OVER(), * "
                "FROM {}".format(self._table)
            )

            if max_lines is None:
                data = cursor.fetchall()
            else:
                data = cursor.fetchmany(max_lines)

            writer.writerows(data)

    def export_sqlite(self, file_name, table_name):
        """
        Export ranking to SQLite database file.
        
        Creates a new SQLite database (or replaces existing table) with
        the ranked results. Useful for further analysis with SQL tools
        or integration with other analysis workflows.
        
        Args:
            file_name (str): Output SQLite database file path
            table_name (str): Name of table to create in database
        """
        self._database.conn.executescript(
            "ATTACH DATABASE '{0}' AS exportdb; "
            "DROP TABLE IF EXISTS exportdb.{1}; "
            "CREATE TABLE 'exportdb'.'{1}' AS "
            "SELECT * "
            "FROM {2}; "
            "DETACH DATABASE exportdb;".format(file_name, table_name, self._table)
        )
