"""
"""
from mpi4py import MPI
import multiprocessing
from enum import IntEnum, auto


class RequestType(IntEnum):
    get_task = 0
    save = auto()
    run_task = auto()
    terminate = auto()


class Message:
    def __init__(self, request, data=None):
        self.request = request
        self.data = data


class MPIPlatform:
    def __init__(self, parameters, pipeline):
        self.local_cpu = parameters["local_cpu"]
        self.pipeline = pipeline
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.root = 0
        self.manager = None
        self.runner = None

    def run(self):
        if self.rank == self.root:
            self.manager = MPIManager(self.comm, self.local_cpu, self.pipeline)
            self.manager.run()

        else:
            self.runner = MPIRunner(self.comm, self.root, self.local_cpu, self.pipeline)
            self.runner.run()

        # Synchronize all runners so we can run again
        self.comm.Barrier()


class MPIManager:
    def __init__(self, mpi_comm, local_cpu, pipeline):
        self.comm = mpi_comm
        self.rank = mpi_comm.Get_rank()
        self.status = MPI.Status()
        self.pipeline = pipeline
        self.alive = mpi_comm.Get_size() - 1 # Number of MPIRunners
        self.not_saved = 0

    def run(self):
        while self.alive or self.not_saved:
            message = self.comm.recv(source=MPI.ANY_SOURCE, status=self.status)

            if message.request == RequestType.get_task:
                self.send_task(message.data)

            if message.request == RequestType.save:
                self.pipeline.save_results(message.data)
                self.not_saved -= len(message.data)

        self.pipeline.finalize()

    def send_task(self, size):
        tasks = self.pipeline.next_job(size)

        if tasks:
            message = Message(RequestType.run_task, tasks)
            self.not_saved += len(tasks)

        else:
            message = Message(RequestType.terminate)
            self.alive -= 1

        runner = self.status.Get_source()
        self.comm.send(message, dest=runner)


class MPIRunner:
    def __init__(self, mpi_comm, root, n_cpu, pipeline):
        self.comm = mpi_comm
        self.rank = mpi_comm.Get_rank()
        self.root = root
        self.n_cpu = n_cpu
        self.pipeline = pipeline

    def run(self):
        if self.n_cpu > 1:
            self.run_multiple()

        else:
            self.run_single()

    def run_multiple(self):
        with multiprocessing.Pool(self.n_cpu, maxtasksperchild=1) as pool:
            running = []
            message = self.request_task(self.n_cpu)

            while message.request != RequestType.terminate:
                for data in message.data:
                    running.append(pool.apply_async(self.pipeline.run_pipeline, (data,)))

                results, ended = self.get_ready(running, 0.1)

                # Un-track finished processes
                for process in ended:
                    running.remove(process)

                self.send_results(results)
                message = self.request_task(len(results))

            results = []
            for process in running:
                results.append(process.get())

            self.send_results(results)

    def get_ready(self, processes, timeout):
        results = []
        ended = []

        while not results:
            for process in processes:
                process.wait(timeout)

                if process.ready():
                    ended.append(process)
                    results.append(process.get())

        return results, ended

    def run_single(self):
        message = self.request_task(1)

        while message.request != RequestType.terminate:
            results = [self.pipeline.run_pipeline(data) for data in message.data]
            self.send_results(results)

            message = self.request_task(1)

    def send_results(self, results):
        self.comm.send(Message(RequestType.save, results), dest=self.root)

    def request_task(self, size):
        return self.comm.sendrecv(
            Message(RequestType.get_task, size),
            dest=self.root,
            source=self.root
        )
