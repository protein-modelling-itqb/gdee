"""
MPI-based distributed execution platform for the GDEE platform.

This module provides a distributed execution platform using MPI (Message Passing Interface)
for parallel processing of protein variants across multiple nodes. It implements a 
manager-worker architecture for efficient load balancing and resource utilization.
"""

from mpi4py import MPI
import multiprocessing
from enum import IntEnum, auto


class RequestType(IntEnum):
    """
    Message types for MPI communication.
    
    Defines the types of messages exchanged between the manager
    and worker nodes during distributed execution.
    """
    get_task = 0      # Worker requests new tasks
    save = auto()     # Worker sends results to save
    run_task = auto() # Manager assigns tasks to worker
    terminate = auto() # Manager signals worker termination


class Message:
    """
    MPI message container for inter-node communication.
    
    Encapsulates request types and associated data for communication
    between manager and worker nodes in the distributed system.
    
    Attributes:
        request (RequestType): Type of request or response
        data: Associated data payload (tasks, results, etc.)
    """
    
    def __init__(self, request, data=None):
        """
        Initialize MPI message.
        
        Args:
            request (RequestType): Message type
            data: Optional data payload
        """
        self.request = request
        self.data = data


class MPIPlatform:
    """
    Distributed execution platform using MPI for GDEE workflows.
    
    This platform coordinates execution across multiple compute nodes using MPI.
    It implements a manager-worker architecture where rank 0 manages job distribution
    and result collection, while other ranks process variants in parallel.
    
    Attributes:
        local_cpu (int): Number of CPU cores per node for local parallelism
        pipeline (Pipeline): Configured pipeline for variant processing
        comm (MPI.Comm): MPI communicator for inter-node communication
        rank (int): MPI rank of current process
        root (int): Root rank (manager node)
        manager (MPIManager): Manager component (rank 0 only)
        runner (MPIRunner): Worker component (non-root ranks only)
    """
    
    def __init__(self, parameters, pipeline):
        """
        Initialize MPI execution platform.
        
        Args:
            parameters (dict): Platform configuration including local CPU count
            pipeline (Pipeline): Configured pipeline for variant processing
        """
        self.local_cpu = parameters["local_cpu"]
        self.pipeline = pipeline
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.root = 0
        self.manager = None
        self.runner = None

    def run(self):
        """
        Execute distributed protein engineering workflow.
        
        Spawns either a manager (rank 0) or worker (other ranks) based on
        MPI rank. The manager coordinates job distribution and result collection,
        while workers process variants through the complete pipeline.
        
        Uses MPI barrier synchronization to ensure all processes complete
        before allowing multiple workflow executions.
        """
        if self.rank == self.root:
            # Root rank acts as job manager
            self.manager = MPIManager(self.comm, self.local_cpu, self.pipeline)
            self.manager.run()

        else:
            # Other ranks act as workers
            self.runner = MPIRunner(self.comm, self.root, self.local_cpu, self.pipeline)
            self.runner.run()

        # Synchronize all processes for potential re-execution
        self.comm.Barrier()


class MPIManager:
    """
    Manager component for distributed MPI execution.
    
    Coordinates job distribution, result collection, and termination signaling
    for the distributed protein engineering workflow. Runs only on rank 0.
    
    Attributes:
        comm (MPI.Comm): MPI communicator
        rank (int): MPI rank (should be 0)
        status (MPI.Status): Status object for message reception
        pipeline (Pipeline): Pipeline for job generation and result saving
        alive (int): Number of active worker nodes
        not_saved (int): Number of jobs sent but not yet saved
    """
    
    def __init__(self, mpi_comm, local_cpu, pipeline):
        """
        Initialize MPI manager.
        
        Args:
            mpi_comm (MPI.Comm): MPI communicator
            local_cpu (int): Local CPU count (unused by manager)
            pipeline (Pipeline): Pipeline for job coordination
        """
        self.comm = mpi_comm
        self.rank = mpi_comm.Get_rank()
        self.status = MPI.Status()
        self.pipeline = pipeline
        self.alive = mpi_comm.Get_size() - 1  # Number of worker nodes
        self.not_saved = 0  # Track pending results

    def run(self):
        """
        Execute manager coordination loop.
        
        Continuously processes requests from workers until all workers
        terminate and all results are saved. Handles job distribution
        and result collection through MPI message passing.
        """
        while self.alive or self.not_saved:
            # Receive request from any worker
            message = self.comm.recv(source=MPI.ANY_SOURCE, status=self.status)

            if message.request == RequestType.get_task:
                # Worker requesting new tasks
                self.send_task(message.data)

            elif message.request == RequestType.save:
                # Worker sending results to save
                self.pipeline.save_results(message.data)
                self.not_saved -= len(message.data)

        # Finalize archiving when all work complete
        self.pipeline.finalize()

    def send_task(self, size):
        """
        Send tasks or termination signal to requesting worker.
        
        Args:
            size (int): Number of tasks requested by worker
        """
        # Get next batch of variants from pipeline
        tasks = self.pipeline.next_job(size)

        if tasks:
            # Send tasks to worker
            message = Message(RequestType.run_task, tasks)
            self.not_saved += len(tasks)

        else:
            # No more tasks - signal termination
            message = Message(RequestType.terminate)
            self.alive -= 1

        # Send response to requesting worker
        runner = self.status.Get_source()
        self.comm.send(message, dest=runner)


class MPIRunner:
    """
    Worker component for distributed MPI execution.
    
    Processes protein variants through the complete pipeline using local
    parallelism (multiprocessing) when multiple CPUs are available.
    Communicates with manager for job requests and result submission.
    
    Attributes:
        comm (MPI.Comm): MPI communicator
        rank (int): MPI rank of this worker
        root (int): Manager rank (0)
        n_cpu (int): Number of local CPU cores for parallel processing
        pipeline (Pipeline): Pipeline for variant processing
    """
    
    def __init__(self, mpi_comm, root, n_cpu, pipeline):
        """
        Initialize MPI worker.
        
        Args:
            mpi_comm (MPI.Comm): MPI communicator
            root (int): Manager rank
            n_cpu (int): Number of local CPU cores
            pipeline (Pipeline): Pipeline for variant processing
        """
        self.comm = mpi_comm
        self.rank = mpi_comm.Get_rank()
        self.root = root
        self.n_cpu = n_cpu
        self.pipeline = pipeline

    def run(self):
        """
        Execute worker processing loop.
        
        Chooses between single-threaded or multi-threaded execution
        based on available CPU cores, then processes variants until
        receiving termination signal from manager.
        """
        if self.n_cpu > 1:
            self.run_multiple()
        else:
            self.run_single()

    def run_multiple(self):
        """
        Execute worker with local multiprocessing parallelism.
        
        Uses multiprocessing.Pool to process multiple variants simultaneously
        on available CPU cores. Manages asynchronous job completion and
        result collection while maintaining communication with manager.
        """
        with multiprocessing.Pool(self.n_cpu, maxtasksperchild=1) as pool:
            running = []  # Track running async jobs
            message = self.request_task(self.n_cpu)

            while message.request != RequestType.terminate:
                # Submit new tasks to process pool
                for data in message.data:
                    async_result = pool.apply_async(self.pipeline.run_pipeline, (data,))
                    running.append(async_result)

                # Check for completed jobs with short timeout
                results, ended = self.get_ready(running, 0.1)

                # Remove completed jobs from tracking
                for process in ended:
                    running.remove(process)

                # Send completed results to manager
                self.send_results(results)
                
                # Request new tasks to replace completed ones
                message = self.request_task(len(results))

            # Process any remaining jobs before termination
            results = []
            for process in running:
                results.append(process.get())

            self.send_results(results)

    def get_ready(self, processes, timeout):
        """
        Check for completed multiprocessing jobs.
        
        Args:
            processes (list): List of AsyncResult objects
            timeout (float): Timeout for job completion checks
            
        Returns:
            tuple: (completed_results, completed_processes)
        """
        results = []
        ended = []

        # Wait for at least one job to complete
        while not results:
            for process in processes:
                process.wait(timeout)

                if process.ready():
                    ended.append(process)
                    results.append(process.get())

        return results, ended

    def run_single(self):
        """
        Execute worker with single-threaded processing.
        
        Processes variants sequentially without local parallelism.
        Suitable for memory-constrained environments or when CPU
        cores are limited.
        """
        message = self.request_task(1)

        while message.request != RequestType.terminate:
            # Process variants sequentially through pipeline
            results = [self.pipeline.run_pipeline(data) for data in message.data]
            
            # Send results to manager
            self.send_results(results)
            
            # Request next batch
            message = self.request_task(1)

    def send_results(self, results):
        """
        Send completed results to manager for database storage.
        
        Args:
            results (list): List of completed job data containers
        """
        self.comm.send(Message(RequestType.save, results), dest=self.root)

    def request_task(self, size):
        """
        Request new tasks from manager.
        
        Args:
            size (int): Number of tasks to request
            
        Returns:
            Message: Response from manager with tasks or termination signal
        """
        return self.comm.sendrecv(
            Message(RequestType.get_task, size),
            dest=self.root,
            source=self.root
        )
