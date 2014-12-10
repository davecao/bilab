
from threading import Thread, Condition, Event
try:
    from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty


__all__ = ['ThreadPool', 'Worker']

class Worker(Thread):
    """
    Thread executing tasks from a given tasks queue
    """
    def __init__(self, tasks, resultsQueue):
        super(Worker, self).__init__()
        self.tasks = tasks
        self.resultsQueue = resultsQueue
        self._dismissed = Event()
        self.daemon = True
        self.start()
    
    def run(self):
        """
        Repeatedly process the job queue until told to exit.
        """
        while not self._dismissed.isSet():
            # Thread blocks here, if queue is empty
            func, args, kargs = self.tasks.get()
            try: 
                func(*args, **kargs)
            except Exception, e: 
                print e
            self.tasks.task_done()

    def dismiss(self):
        """
        Set a flag to tell the thread to exit when done with current job.
        """
        self._dismissed.set()

class ThreadPool(object):
    """
    Pool of threads consuming tasks from a queue
    """
    
    def __init__(self, num_threads):
        super(ThreadPool, self).__init__()
        self.tasks = Queue(num_threads)
        self.resultsQueue = Queue()
        self.workers = []

        for _ in range(num_threads): 
            self.workers.append(Worker(self.tasks, self.resultsQueue))

    def dismissWokers(self, num_threads):
        """
        Tell num_threads worker threads to to quit when they're done.
        """
        for i in range(min(num_threads, len(self.workers))):
            worker = self.workers.pop()
            worker.dismiss()

    def add_task(self, func, *args, **kargs):
        """
        Add a task to the queue
        """
        self.tasks.put((func, args, kargs))

    def wait_completion(self):
        """
        Wait for completion of all the tasks in the queue
        """
        self.tasks.join()
