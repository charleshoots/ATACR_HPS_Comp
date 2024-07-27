def Fork(fn, *args, **kwargs):
  """Submits a job to the concurrent.futures ThreadPoolExecutor.

  Args:
    fn: The function to be executed.
    *args: The arguments to be passed to the function.
    **kwargs: The keyword arguments to be passed to the function.

  Returns:
    A Future object representing the job.
  """
  max_workers = kwargs['max_workers']
  with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
    future = executor.submit(fn, *args, **kwargs)
    print('Waiting for tasks to complete')
    wait([future])
    return future