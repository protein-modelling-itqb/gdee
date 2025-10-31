"""
Single-machine execution platform for the GDEE platform.

This module provides a simple, sequential execution platform that processes
protein variants one at a time on a single machine. It's suitable for
small-scale studies and development work.
"""


class SimplePlatform:
    """
    Sequential execution platform for the GDEE platform.

    This platform executes the complete GDEE pipeline sequentially on a single
    machine. It processes variants one by one through the entire workflow:
    variant generation → modeling → quality assessment → docking → measurements.
    
    Attributes:
        pipeline (Pipeline): Configured pipeline with all processing components
    """
    
    def __init__(self, parameters, pipeline):
        """
        Initialize simple execution platform.
        
        Args:
            parameters (dict): Platform configuration (unused for simple platform)
            pipeline (Pipeline): Configured pipeline for variant processing
        """
        self.pipeline = pipeline

    def run(self):
        """
        Execute the protein engineering workflow sequentially.
        
        This method:
        1. Requests variant batches from the pipeline
        2. Processes each variant through all pipeline stages
        3. Saves results to database and archives successful jobs
        4. Continues until no more variants are available
        5. Finalizes archiving and cleanup operations
        
        The platform processes variants in single batches to minimize
        memory usage while maintaining simplicity.
        """
        job_data = self.pipeline.next_job(1)

        while job_data:
            # Process variants through complete pipeline
            results = [self.pipeline.run_pipeline(data) for data in job_data]
            
            # Save results to database and archive files
            self.pipeline.save_results(results)
            
            # Get next batch of variants
            job_data = self.pipeline.next_job(1)

        # Complete archiving and cleanup
        self.pipeline.finalize()
