"""
File archiving system for GDEE platform results.

This module provides archiving functionality for storing and organizing results
from GDEE platform. It creates compressed tar archives in a
buffered, asynchronous manner to optimize I/O performance during pipeline execution.
"""

import tarfile
from path import Path
import multiprocessing as mp


class Archiver:
    """
    Asynchronous file archiver for GDEE platform results.

    This class manages the archiving of job directories and results into compressed
    tar files. It uses a buffering system with background processes to minimize
    I/O blocking during pipeline execution. Archives are created incrementally
    based on the buffer size to manage disk space and file handling.
    
    The archiver is used by the Pipeline class to store successful job results.
    
    
    Attributes:
        archive (Path): Base name for archive files
        format (str): Format string for archive file naming (includes .tar extension)
        buffer_size (int): Maximum number of files to buffer before flushing
        _buffer (list): Internal buffer storing (path, new_name) tuples
        _job (Process): Background multiprocessing job for tar creation
        _index (int): Current archive file index number
    """
    
    def __init__(self, name, format, buffer_size):
        """
        Initialize the file archiver.
        
        Finds the next available archive index by checking existing files
        to avoid overwriting previous archives.
        
        Args:
            name (str): Base name for archive files (e.g., "files")
            format (str): Format string for numbering (e.g., ".{:06d}")
            buffer_size (int): Number of files to buffer before creating archive
        """
        self.archive = Path(name)
        self.format = format + ".tar"
        self.buffer_size = buffer_size
        self._buffer = []
        self._job = None

        # Find next available archive index
        self._index = 0
        while self._next_filename().exists():
            pass
        self._index -= 1

    def _next_filename(self):
        """
        Generate the next archive filename.
        
        Creates filenames like "files.000001.tar", "files.000002.tar", etc.
        based on the format string provided during initialization.
        
        Returns:
            Path: Next archive file path
        """
        new_file = self.archive + self.format.format(self._index)
        self._index += 1
        return new_file

    def add(self, path, new_name):
        """
        Add a file or directory to the archive buffer.
        
        If the buffer reaches its maximum size, automatically triggers
        a flush operation to create the current archive and start a new one.
        
        Args:
            path (str or Path): Source path to archive
            new_name (str): Name to use within the archive
        """
        if len(self._buffer) >= self.buffer_size:
            self.flush()

        self._buffer.append((path, new_name))

    def finalize(self):
        """
        Complete all archiving operations.
        
        Flushes any remaining buffered files and waits for background
        processes to complete. This should be called at the end of
        pipeline execution to ensure all results are properly archived.
        """
        self.flush()
        self.join()

    def join(self):
        """
        Wait for the current background archiving process to complete.
        
        This method blocks until any running multiprocessing job finishes,
        ensuring that archive creation is complete before proceeding.
        """
        if self._job is not None:
            self._job.join()

    def flush(self):
        """
        Create an archive from the current buffer contents.
        
        Waits for any existing background job to complete, then starts
        a new multiprocessing job to create a tar archive containing
        all buffered files. The buffer is cleared after starting the job.
        
        If the buffer is empty, no action is taken.
        """
        self.join()

        if not self._buffer:
            return

        self._job = mp.Process(
            target=self._save_results, 
            args=(self._next_filename(), self._buffer)
        )
        self._job.start()
        self._buffer = []

    def _save_results(self, file_name, data):
        """
        Background process worker for creating tar archives.
        
        This method runs in a separate process to avoid blocking the main
        pipeline execution. It creates a tar archive and adds all files
        from the data list, then removes the source directories to save
        disk space.
        
        Args:
            file_name (Path): Target archive file path
            data (list): List of (source_path, archive_name) tuples to archive
        """
        with tarfile.open(file_name, "a") as tar:
            while data:
                path, new_name = data.pop()
                tar.add(path, new_name)
                Path(path).rmtree_p()
