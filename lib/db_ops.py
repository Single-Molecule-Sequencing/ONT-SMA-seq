"""Database operations for CentralDB and SMADB."""

import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Any

# Column whitelists for SQL injection prevention
EXPERIMENTS_COLUMNS = {
    'experiment_path', 'instrument', 'flow_cell_id', 'sample_id',
    'protocol', 'started', 'pod5_count', 'pod5_dir', 'created_at'
}

BASECALL_RUNS_COLUMNS = {
    'model_tier', 'model_version', 'trim', 'mod_bitflag', 'dorado_version',
    'dorado_args', 'batch_size', 'emit_moves', 'gpu_model', 'slurm_job_id',
    'runtime_seconds', 'bam_path', 'sma_db_path', 'status', 'created_at'
}

LIBRARY_SPECS_COLUMNS = {
    'lib_version', 'description', 'created_by', 'project', 'pore_type',
    'flow_cell_type', 'sequencing_kit', 'status', 'created_at'
}

READS_COLUMNS = {
    'readseq', 'readlen', 'q_bc', 'ed', 'q_ld', 'ref_id', 'end_reason'
}

RUN_SUMMARIES_COLUMNS = {
    'total_reads', 'matched_reads', 'mean_qbc', 'median_qbc', 'mean_ed',
    'median_ed', 'mean_qld', 'er_signal_pos', 'er_signal_neg', 'er_unblock',
    'er_other', 'created_at'
}


def _validate_columns(kwargs_keys: set, valid_columns: set, table_name: str) -> None:
    """Validate that all kwargs keys are valid column names.

    Args:
        kwargs_keys: Set of column names from kwargs
        valid_columns: Set of valid column names for the table
        table_name: Name of the table (for error messages)

    Raises:
        ValueError: If any invalid column names are found
    """
    invalid_cols = kwargs_keys - valid_columns
    if invalid_cols:
        raise ValueError(f"Invalid columns for {table_name} table: {invalid_cols}")


class CentralDB:
    """Operations for central nanopore_unified.db database.

    Provides methods for managing experiments, basecall runs, library specs,
    and their relationships.

    Example usage:
        with CentralDB('nanopore_unified.db') as db:
            db.insert_experiment('exp1', sample_id='sample1')
            exp = db.get_experiment('exp1')
    """

    def __init__(self, db_path: str):
        """Initialize connection to central database.

        Args:
            db_path: Path to the central database file
        """
        self.db_path = db_path
        self._conn: Optional[sqlite3.Connection] = None
        self._connect()

    def _connect(self) -> None:
        """Establish database connection with row factory."""
        self._conn = sqlite3.connect(self.db_path)
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA foreign_keys = ON;")

    @property
    def conn(self) -> sqlite3.Connection:
        """Get the active database connection."""
        if self._conn is None:
            raise RuntimeError("Database connection is closed")
        return self._conn

    def close(self) -> None:
        """Close the database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self) -> 'CentralDB':
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit - rollback on exception and close connection."""
        if exc_type is not None and self._conn is not None:
            self._conn.rollback()
        self.close()

    def _row_to_dict(self, row: Optional[sqlite3.Row]) -> Optional[Dict[str, Any]]:
        """Convert sqlite3.Row to dict."""
        if row is None:
            return None
        return dict(row)

    def _rows_to_list(self, rows: List[sqlite3.Row]) -> List[Dict[str, Any]]:
        """Convert list of sqlite3.Row to list of dicts."""
        return [dict(row) for row in rows]

    # ========== Experiment Operations ==========

    def insert_experiment(self, exp_id: str, **kwargs) -> None:
        """Insert a new experiment record.

        Args:
            exp_id: Unique experiment identifier
            **kwargs: Optional fields matching experiments table schema:
                - experiment_path: Path to experiment directory
                - instrument: Instrument type (e.g., 'P2S', 'MinION')
                - flow_cell_id: Flow cell identifier
                - sample_id: Sample identifier
                - protocol: Protocol name
                - started: Start timestamp
                - pod5_count: Number of POD5 files
                - pod5_dir: Path to POD5 directory

        Raises:
            ValueError: If invalid column names are provided
            sqlite3.IntegrityError: If exp_id already exists
        """
        # Validate column names to prevent SQL injection
        _validate_columns(set(kwargs.keys()), EXPERIMENTS_COLUMNS, 'experiments')

        columns = ['exp_id'] + list(kwargs.keys())
        placeholders = ['?'] * len(columns)
        values = [exp_id] + list(kwargs.values())

        sql = f"INSERT INTO experiments ({', '.join(columns)}) VALUES ({', '.join(placeholders)})"
        try:
            self.conn.execute(sql, values)
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def get_experiment(self, exp_id: str) -> Optional[Dict[str, Any]]:
        """Retrieve an experiment by ID.

        Args:
            exp_id: Experiment identifier

        Returns:
            Dict with experiment fields or None if not found
        """
        cursor = self.conn.execute(
            "SELECT * FROM experiments WHERE exp_id = ?",
            (exp_id,)
        )
        return self._row_to_dict(cursor.fetchone())

    def list_experiments(self) -> List[Dict[str, Any]]:
        """List all experiments.

        Returns:
            List of dicts with experiment fields
        """
        cursor = self.conn.execute(
            "SELECT * FROM experiments ORDER BY created_at DESC"
        )
        return self._rows_to_list(cursor.fetchall())

    # ========== Basecall Run Operations ==========

    def insert_basecall_run(self, run_id: str, exp_id: str, **kwargs) -> None:
        """Insert a new basecall run record.

        Args:
            run_id: Unique run identifier
            exp_id: Parent experiment identifier
            **kwargs: Optional fields matching basecall_runs table schema:
                - model_tier: 'fast', 'hac', or 'sup'
                - model_version: Model version string
                - trim: Trim setting (0/1)
                - mod_bitflag: Modification bitflag
                - dorado_version: Dorado version
                - dorado_args: Full dorado command arguments
                - batch_size: Batch size used
                - emit_moves: Whether moves were emitted (0/1)
                - gpu_model: GPU model used
                - slurm_job_id: SLURM job ID
                - runtime_seconds: Total runtime
                - bam_path: Output BAM path
                - sma_db_path: Path to SMA database
                - status: Run status ('pending', 'running', 'complete', 'failed')

        Raises:
            ValueError: If invalid column names are provided
            sqlite3.IntegrityError: If run_id exists or exp_id not found
        """
        # Validate column names to prevent SQL injection
        _validate_columns(set(kwargs.keys()), BASECALL_RUNS_COLUMNS, 'basecall_runs')

        columns = ['run_id', 'exp_id'] + list(kwargs.keys())
        placeholders = ['?'] * len(columns)
        values = [run_id, exp_id] + list(kwargs.values())

        sql = f"INSERT INTO basecall_runs ({', '.join(columns)}) VALUES ({', '.join(placeholders)})"
        try:
            self.conn.execute(sql, values)
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def get_basecall_run(self, run_id: str) -> Optional[Dict[str, Any]]:
        """Retrieve a basecall run by ID.

        Args:
            run_id: Run identifier

        Returns:
            Dict with run fields or None if not found
        """
        cursor = self.conn.execute(
            "SELECT * FROM basecall_runs WHERE run_id = ?",
            (run_id,)
        )
        return self._row_to_dict(cursor.fetchone())

    def list_basecall_runs(
        self,
        exp_id: Optional[str] = None,
        status: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """List basecall runs with optional filtering.

        Args:
            exp_id: Filter by experiment ID (optional)
            status: Filter by status (optional)

        Returns:
            List of dicts with run fields
        """
        conditions = []
        params = []

        if exp_id is not None:
            conditions.append("exp_id = ?")
            params.append(exp_id)

        if status is not None:
            conditions.append("status = ?")
            params.append(status)

        sql = "SELECT * FROM basecall_runs"
        if conditions:
            sql += " WHERE " + " AND ".join(conditions)
        sql += " ORDER BY created_at DESC"

        cursor = self.conn.execute(sql, params)
        return self._rows_to_list(cursor.fetchall())

    def update_basecall_run_status(
        self,
        run_id: str,
        status: str,
        sma_db_path: Optional[str] = None
    ) -> None:
        """Update basecall run status and optionally SMA database path.

        Args:
            run_id: Run identifier
            status: New status value
            sma_db_path: Optional path to SMA database
        """
        try:
            if sma_db_path is not None:
                self.conn.execute(
                    "UPDATE basecall_runs SET status = ?, sma_db_path = ? WHERE run_id = ?",
                    (status, sma_db_path, run_id)
                )
            else:
                self.conn.execute(
                    "UPDATE basecall_runs SET status = ? WHERE run_id = ?",
                    (status, run_id)
                )
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    # ========== Library Spec Operations ==========

    def insert_library_spec(self, lib_id: str, lib_name: str, **kwargs) -> None:
        """Insert a new library specification.

        Args:
            lib_id: Unique library identifier
            lib_name: Library name (required)
            **kwargs: Optional fields matching library_specs table schema:
                - lib_version: Version string
                - description: Library description
                - created_by: Creator name
                - project: Project name
                - pore_type: Pore type (e.g., 'R10.4.1')
                - flow_cell_type: Flow cell type
                - sequencing_kit: Kit used
                - status: 'draft', 'active', 'archived', 'deprecated'

        Raises:
            ValueError: If invalid column names are provided
            sqlite3.IntegrityError: If lib_id already exists
        """
        # Validate column names to prevent SQL injection
        _validate_columns(set(kwargs.keys()), LIBRARY_SPECS_COLUMNS, 'library_specs')

        columns = ['lib_id', 'lib_name'] + list(kwargs.keys())
        placeholders = ['?'] * len(columns)
        values = [lib_id, lib_name] + list(kwargs.values())

        sql = f"INSERT INTO library_specs ({', '.join(columns)}) VALUES ({', '.join(placeholders)})"
        try:
            self.conn.execute(sql, values)
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def get_library_spec(self, lib_id: str) -> Optional[Dict[str, Any]]:
        """Retrieve a library specification by ID.

        Args:
            lib_id: Library identifier

        Returns:
            Dict with library spec fields or None if not found
        """
        cursor = self.conn.execute(
            "SELECT * FROM library_specs WHERE lib_id = ?",
            (lib_id,)
        )
        return self._row_to_dict(cursor.fetchone())

    def link_library_experiment(
        self,
        lib_id: str,
        exp_id: str,
        relationship: str = 'primary',
        run_purpose: Optional[str] = None
    ) -> None:
        """Link a library specification to an experiment.

        Args:
            lib_id: Library identifier
            exp_id: Experiment identifier
            relationship: Type of relationship (default: 'primary')
            run_purpose: Purpose of the run (optional)

        Raises:
            sqlite3.IntegrityError: If link already exists or IDs not found
        """
        try:
            self.conn.execute(
                """INSERT INTO lib_experiment_link (lib_id, exp_id, relationship, run_purpose)
                   VALUES (?, ?, ?, ?)""",
                (lib_id, exp_id, relationship, run_purpose)
            )
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    # ========== Query Operations ==========

    def get_pending_sma_analysis(self) -> List[Dict[str, Any]]:
        """Get basecall runs ready for SMA analysis.

        Returns runs where status='complete' and sma_db_path is NULL.

        Returns:
            List of dicts with run fields
        """
        cursor = self.conn.execute(
            """SELECT * FROM basecall_runs
               WHERE status = 'complete' AND sma_db_path IS NULL
               ORDER BY created_at"""
        )
        return self._rows_to_list(cursor.fetchall())


class SMADB:
    """Operations for per-experiment SMA_{exp_id}.db databases.

    Provides methods for storing and retrieving per-read metrics
    with performance optimizations for bulk inserts.

    Example usage:
        with SMADB('SMA_exp1.db') as db:
            db.insert_read('read1', 'run1', readlen=5000, q_bc=25.5)
            summary = db.compute_run_summary('run1')
    """

    def __init__(self, db_path: str):
        """Initialize connection to SMA database.

        Args:
            db_path: Path to the SMA database file

        Note:
            Applies performance tuning PRAGMAs for bulk operations.
        """
        self.db_path = db_path
        self._conn: Optional[sqlite3.Connection] = None
        self._connect()

    def _connect(self) -> None:
        """Establish database connection with performance tuning."""
        self._conn = sqlite3.connect(self.db_path)
        self._conn.row_factory = sqlite3.Row
        # Performance tuning for bulk inserts
        self._conn.execute("PRAGMA synchronous = OFF;")
        self._conn.execute("PRAGMA journal_mode = MEMORY;")

    @property
    def conn(self) -> sqlite3.Connection:
        """Get the active database connection."""
        if self._conn is None:
            raise RuntimeError("Database connection is closed")
        return self._conn

    def close(self) -> None:
        """Close the database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self) -> 'SMADB':
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit - rollback on exception and close connection."""
        if exc_type is not None and self._conn is not None:
            self._conn.rollback()
        self.close()

    def _row_to_dict(self, row: Optional[sqlite3.Row]) -> Optional[Dict[str, Any]]:
        """Convert sqlite3.Row to dict."""
        if row is None:
            return None
        return dict(row)

    # ========== Read Operations ==========

    def insert_read(self, read_id: str, run_id: str, **kwargs) -> None:
        """Insert a single read record.

        Args:
            read_id: Unique read identifier
            run_id: Parent run identifier
            **kwargs: Optional fields matching reads table schema:
                - readseq: Read sequence
                - readlen: Read length
                - q_bc: Quality score from basecalling
                - ed: Edit distance to reference
                - q_ld: Quality score from length distribution
                - ref_id: Matched reference ID
                - end_reason: End reason from POD5

        Raises:
            ValueError: If invalid column names are provided
            sqlite3.IntegrityError: If (read_id, run_id) already exists
        """
        # Validate column names to prevent SQL injection
        _validate_columns(set(kwargs.keys()), READS_COLUMNS, 'reads')

        columns = ['read_id', 'run_id'] + list(kwargs.keys())
        placeholders = ['?'] * len(columns)
        values = [read_id, run_id] + list(kwargs.values())

        sql = f"INSERT INTO reads ({', '.join(columns)}) VALUES ({', '.join(placeholders)})"
        try:
            self.conn.execute(sql, values)
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def insert_reads_batch(self, reads: List[tuple]) -> None:
        """Bulk insert read records for performance.

        Args:
            reads: List of tuples, each containing:
                (read_id, run_id, readseq, readlen, q_bc, ed, q_ld, ref_id, end_reason)
                Use None for optional fields.

        Note:
            This method is optimized for bulk inserts and commits
            only once at the end of all inserts. Column validation is not
            performed since this uses a fixed schema.
        """
        sql = """INSERT INTO reads
                 (read_id, run_id, readseq, readlen, q_bc, ed, q_ld, ref_id, end_reason)
                 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)"""
        try:
            self.conn.executemany(sql, reads)
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    # ========== Run Summary Operations ==========

    def insert_run_summary(self, run_id: str, **kwargs) -> None:
        """Insert or replace a run summary record.

        Args:
            run_id: Run identifier
            **kwargs: Optional fields matching run_summaries table schema:
                - total_reads: Total read count
                - matched_reads: Reads matched to references
                - mean_qbc: Mean quality score
                - median_qbc: Median quality score
                - mean_ed: Mean edit distance
                - median_ed: Median edit distance
                - mean_qld: Mean quality from length distribution
                - er_signal_pos: Count of signal_positive end reasons
                - er_signal_neg: Count of signal_negative end reasons
                - er_unblock: Count of unblock end reasons
                - er_other: Count of other end reasons

        Raises:
            ValueError: If invalid column names are provided
        """
        # Validate column names to prevent SQL injection
        _validate_columns(set(kwargs.keys()), RUN_SUMMARIES_COLUMNS, 'run_summaries')

        columns = ['run_id'] + list(kwargs.keys())
        placeholders = ['?'] * len(columns)
        values = [run_id] + list(kwargs.values())

        sql = f"""INSERT OR REPLACE INTO run_summaries
                  ({', '.join(columns)}) VALUES ({', '.join(placeholders)})"""
        try:
            self.conn.execute(sql, values)
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def get_run_summary(self, run_id: str) -> Optional[Dict[str, Any]]:
        """Retrieve a run summary by ID.

        Args:
            run_id: Run identifier

        Returns:
            Dict with summary fields or None if not found
        """
        cursor = self.conn.execute(
            "SELECT * FROM run_summaries WHERE run_id = ?",
            (run_id,)
        )
        return self._row_to_dict(cursor.fetchone())

    def compute_run_summary(self, run_id: str) -> Dict[str, Any]:
        """Compute summary statistics from reads table.

        Aggregates metrics for all reads in a run and returns
        computed statistics without storing them.

        Args:
            run_id: Run identifier

        Returns:
            Dict with computed summary statistics:
                - total_reads: Total read count
                - matched_reads: Reads with ref_id set
                - mean_qbc: Mean quality score
                - median_qbc: Median quality score (approximated)
                - mean_ed: Mean edit distance
                - median_ed: Median edit distance (approximated)
                - mean_qld: Mean quality from length distribution
                - er_signal_pos: signal_positive count
                - er_signal_neg: signal_negative count
                - er_unblock: unblock_mux_change count
                - er_other: Other end reason count
        """
        # Get aggregate statistics
        cursor = self.conn.execute("""
            SELECT
                COUNT(*) as total_reads,
                SUM(CASE WHEN ref_id IS NOT NULL THEN 1 ELSE 0 END) as matched_reads,
                AVG(q_bc) as mean_qbc,
                AVG(ed) as mean_ed,
                AVG(q_ld) as mean_qld,
                SUM(CASE WHEN end_reason = 'signal_positive' THEN 1 ELSE 0 END) as er_signal_pos,
                SUM(CASE WHEN end_reason = 'signal_negative' THEN 1 ELSE 0 END) as er_signal_neg,
                SUM(CASE WHEN end_reason = 'unblock_mux_change' THEN 1 ELSE 0 END) as er_unblock,
                SUM(CASE WHEN end_reason NOT IN ('signal_positive', 'signal_negative', 'unblock_mux_change')
                         AND end_reason IS NOT NULL THEN 1 ELSE 0 END) as er_other
            FROM reads
            WHERE run_id = ?
        """, (run_id,))

        row = cursor.fetchone()
        result = dict(row) if row else {}

        # Compute median q_bc (SQLite doesn't have native MEDIAN)
        if result.get('total_reads', 0) > 0:
            cursor = self.conn.execute("""
                SELECT q_bc FROM reads
                WHERE run_id = ? AND q_bc IS NOT NULL
                ORDER BY q_bc
                LIMIT 1 OFFSET (
                    SELECT COUNT(*) / 2 FROM reads
                    WHERE run_id = ? AND q_bc IS NOT NULL
                )
            """, (run_id, run_id))
            median_row = cursor.fetchone()
            result['median_qbc'] = median_row['q_bc'] if median_row else None

            # Compute median ed
            cursor = self.conn.execute("""
                SELECT ed FROM reads
                WHERE run_id = ? AND ed IS NOT NULL
                ORDER BY ed
                LIMIT 1 OFFSET (
                    SELECT COUNT(*) / 2 FROM reads
                    WHERE run_id = ? AND ed IS NOT NULL
                )
            """, (run_id, run_id))
            median_row = cursor.fetchone()
            result['median_ed'] = median_row['ed'] if median_row else None
        else:
            result['median_qbc'] = None
            result['median_ed'] = None

        return result
