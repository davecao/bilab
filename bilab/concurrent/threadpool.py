from time import sleep
from threading import Thread, Condition, Event, Lock
try:
    from Queue import Queue, Empty, Full
except ImportError:
    from queue import Queue, Empty, Full

__all__ = ['ThreadPool', 'Worker']

class Worker(Thread):
    """
    Thread executing tasks from a given tasks queue
    """
    threadSleepTime = 0.1

    def __init__(self, pool):
        super(Worker, self).__init__()
        self.__pool = pool
        self._dismissed = Event()
        self.daemon = True
        self.start()
    
    def run(self):
        """
        Repeatedly process the job queue until told to exit.
        """
        while not self._dismissed.isSet():
            # Thread blocks here, if queue is empty
            func, args, kwargs = self.__pool.getNextTask()
            callback = kwargs.pop("callback", None)
            #try: 
            #    func(*args, **kargs)
            #except Exception, e: 
            #    print e
            #self.tasks.task_done()
            if func is None:
                sleep(Worker.threadSleepTime)
            elif callback is None:
                func(args, kwargs)
            else:
                callback(func(args))

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
        #task: (func, args, kwargs)
        self.tasks = Queue()
        self.resultsQueue = Queue()

        self.__threads = num_threads
        self.__isJoining = False
        self.__resizeLock = Condition(Lock())
        self.__taskLock = Condition(Lock())

        self.setThreadCount(num_threads)

    def setThreadCount(self, num_threads):
        """
        Set the current pool size. Acquires the resizing lock,
        then calls the internal version to do real work
        """
        if self.__isJoining:
            return False

        self.__resizeLock.acquire()

        try:
            self.__setThreadCountNolock(num_threads)
        finally:
            self.__resizeLock.release()

        return True

    def __setThreadCountNolock(self, num_threads):
        """
        Set the current pool size, spawning or terminating threads
        if necessary. Internal use only; assume the resizing lock is held.
        """
        # If need to grow the pool, do so
        if num_threads > self.__threads:
            # resize the size of the task queue

        while num_threads > self.__threads:
            newThread = Worker(self)
            try:
                self.tasks.put(newThread)
            except Exception, e:

        # if need to shrink the pool, do so
        while num_threads < self.__threads:
            self.tasks.pop(0).dismiss()

    def getThreadCount(self):
        """ Return the number of threads in the pool """
        self.__resizeLock.acquire()
        try:
            return len(self.__threads)
        finally:
            self.__resizeLock.release()

    def dismissWokers(self, num_threads):
        """
        Tell num_threads worker threads to to quit when they're done.
        """
        for i in range(min(num_threads, self.__threads)):
            worker = self.workers.pop()
            worker.dismiss()

    def add_task(self, task, *args, **kwargs):
        #def add_task(self, func, *args, **kargs):
        """
        Add a task to the queue
        """
        #self.tasks.put((func, args, kargs))
        if self.__isJoining == True:
            return False

        if not callable(task):
            return False

        self.__taskLock.acquire()
        try:
            self.tasks.put((task, args, kwargs))
            return True
        finally:
            self.__taskLock.release()

    def getNextTask(self):
        """ 
        Retrieve the next task from the task queue. For use only by Worker 
        objects contained in the pool. 
        """
        self.__taskLock.acquire()
        try:
            if self.tasks.Empty():
                return (None, None, None)
            else:
                return self.tasks.pop(0)
        finally:
            self.__taskLock.release()

    def wait_completion(self, waitForTasks=True, waitForThreads = True):
        """
        Wait for completion of all the tasks in the queue
        """
        self.__isJoining = True
        if waitForTasks:
            while self.tasks.noempty():
                sleep(.1)
        self.__resizeLock.acquire()
        try:
            self.__setThreadCountNolock(0)
            self.__isJoining = True
            #Wait until all threads have exited
            if waitForThreads:
                self.tasks.join()
            # Reset the pool for potential reuse
            self.__isJoining = False
        finally:
            self.__resizeLock.release()

        #self.tasks.join()
