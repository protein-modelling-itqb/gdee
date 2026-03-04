"""Database filtering and ranking for result analysis."""


import csv
import uuid
import numbers
from ..database import Database


def _generate_table_name():
    """Generate unique temporary table name.
    
    Returns:
        str: Unique table name
    """
    return "R" + str(uuid.uuid4()).replace("-", "")


class Metric:
    """Represents a measurement metric in the database."""
    def __init__(self, name, database):
        """Initialize metric reference.
        
        Args:
            name: Metric name
            database: Database instance
        """
        self._name = name
        self._database = database
        self._metric_id = database.fetch_metric_id(name)
        if self._metric_id is None:
            raise RuntimeError("Metric '{}' not found".format(name))

    @property
    def name(self):
        return self._name

    def __lt__(self, value):
        return self._compare("<", value)

    def __le__(self, value):
        return self._compare("<=", value)

    def __gt__(self, value):
        return self._compare(">", value)

    def __ge__(self, value):
        return self._compare(">=", value)

    def __eq__(self, value):
        return self._compare("=", value)

    def __ne__(self, value):
        return self._compare("!=", value)

    def _compare(self, operator, value):
        if isinstance(value, numbers.Number):
            return self._compare_number(operator, value)

        elif isinstance(value, Metric):
            return self._compare_metrics(operator, value)

    def _compare_number(self, operator, value):
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
    """Represents a filter rule for pose selection."""
    def __init__(self, table, database):
        """Initialize filter rule.
        
        Args:
            table: Temporary table name with pose_ids
            database: Database instance
        """
        self._table = table
        self._database = database

    def __and__(self, other):
        """Combine rules with AND operation.
        
        Args:
            other: Another Rule object
            
        Returns:
            Rule: New rule with intersection
        """
        table = _generate_table_name()
        self._database.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT pose_id "
            "FROM {} "
            "WHERE pose_id IN (SELECT pose_id FROM {});".format(table, self._table, other._table))
        return Rule(table, self._database)

    def __or__(self, other):
        """Combine rules with OR operation.
        
        Args:
            other: Another Rule object
            
        Returns:
            Rule: New rule with union
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
        """Invert rule with NOT operation.
        
        Returns:
            Rule: Inverted rule
        """
        table = _generate_table_name()
        self._database.conn.execute(
            "CREATE TEMP TABLE {} AS "
            "SELECT pose_id "
            "FROM Measurements "
            "WHERE pose_id NOT IN (SELECT pose_id FROM {});".format(table, self._table))
        return Rule(table, self._database)

    def __bool__(self):
        """Prevent boolean conversion.
        
        Raises:
            TypeError: Rule cannot be converted to bool
        """
        raise TypeError("Rule not convertible to bool. Use the bitwise operators ~ (not), & (and), | (or) for boolean expressions.")

    def rank(self, ascending=True):
        return Rank.from_rule(self._table, self._database, ascending)


class Rank:
    """Represents ranked results from database query."""
    def __init__(self, table, database):
        """Initialize ranked results.
        
        Args:
            table: Table name with ranked poses
            database: Database instance
        """
        self._table = table
        self._database = database

    @staticmethod
    def from_rule(poses_table, database, ascending):
        """Create ranked results from a rule.
        
        Args:
            poses_table: Table with selected pose_ids
            database: Database instance
            ascending: Sort ascending if True
            
        Returns:
            Rank: Ranked results object
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
        """Filter by number of mutations.
        
        Args:
            num_mutations: Number of mutations to filter by
            
        Returns:
            Rank: Filtered ranked results
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
        """Filter to only wildtype variants.
        
        Returns:
            Rank: Filtered ranked results
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
        """Export results to CSV file.
        
        Args:
            file_name: Output CSV file path
            max_lines: Maximum number of rows to export
        """
        with open(file_name, "w") as fd:
            writer = csv.writer(fd)
            writer.writerow(("rank", "energy", "name", "directory", "is_wildtype", "model_id", "variant_id", "eval_id", "pose_index", "docking_file"))
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
        """Export results to SQLite file.
        
        Args:
            file_name: Output database file path
            table_name: Table name in output database
            
        Returns:
            Rank: Results in output database
        """
        self._database.conn.executescript(
            "ATTACH DATABASE '{0}' AS exportdb; "
            "DROP TABLE IF EXISTS exportdb.{1}; "
            "CREATE TABLE 'exportdb'.'{1}' AS "
            "SELECT * "
            "FROM {2}; "
            "DETACH DATABASE exportdb;".format(file_name, table_name, self._table)
        )
