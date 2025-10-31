"""
"""


class SimplePlatform:
    def __init__(self, parameters, pipeline):
        self.pipeline = pipeline

    def run(self):
        job_data = self.pipeline.next_job(1)

        while job_data:
            results = [self.pipeline.run_pipeline(data) for data in job_data]
            self.pipeline.save_results(results)
            job_data = self.pipeline.next_job(1)

        self.pipeline.finalize()
