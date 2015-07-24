# -*- coding: utf-8 -*-

from functools import  wraps
import multiprocessing
import logging
from time import sleep
from threading import Thread, Condition, Event, Lock
try:
    from Queue import Queue, Empty, Full
except ImportError:
    from queue import Queue, Empty, Full

__all__ = ['AbstractThreadPool', 'AbstractLocalThreadPool']

# decorator borrowed from Mozilla mxr
def abstractmethod(method):
    line = method.func_code.co_firstlineno
    filename = method.func_code.co_filename
    @wraps(method)
    def not_implemented(*args, **kwargs):
        raise NotImplementedError('Abstract method %s at File "%s", line %s'
            'should be implemented by a concrete class' %
            (repr(method), filename, line))
    return not_implemented

class AbstractThreadPool(object):
    """ docstring for class

    .. note::

    """

    def __init__(self, *args, **kwargs):
        """Initialization.

        Args:
           *args:

        Kwargs:
            **kwargs:
        """
        super(AbstractThreadPool, self).__init__()
        self.logger = logging.getLogger(__name__)

    @abstractmethod
    def add_task(self, taskid, task, *args, **kwargs):
        """
        Add task into queue
        """

    @abstractmethod
    def wait_completion(self, waitForTasks=True, waitForThreads = True):
        """
        Wait for completion of all the tasks in the queue
        """

class AbstractLocalThreadPool(AbstractThreadPool):
    """
    Base class for pools that rely on local threads and 
    a queue to dispatch jobs
    """
    #How often workers check for new jobs in the queue
    WorkerCheckInterval = 0.5 
    
    def __init__(self, *args, **kwargs):
        """Initialization

        Args:
            num_threads(int): 
                number of threads, defaults to number of cores 
                that are available on system
        """
        super(AbstractLocalThreadPool, self).__init__(*args, **kwargs)
        
        # Store tasks
        self.tasks = Queue()
        #Threads put themselves here when they die
        self.dead_thread_queue = Queue()
        
        self.shutdown_event = Event()
        self.shutdown_event.clear()

        self.WorkerCheckInterval=kwargs.get('WorkerCheckInterval', 0.5)
        self._max_threads = kwargs.get('num_threads', 
                                        multiprocessing.cpu_count())
        if self._max_threads is None:
            self._max_threads = multiprocessing.cpu_count()
        
        self.keep_alive_thread = None
        self.threads = []

    @abstractmethod
    def add_worker_thread(self):
        """
        Implemented by derived class and return a thread object
        """

    @abstractmethod
    def add_task(self, taskid, task, *args, **kwargs):
        """
        Add task into queue
        """

    @abstractmethod
    def shutdown(self):
        """ Shutdown the pool """

    def __start_keep_alive_thread(self):
        """ Keep alive threads """
        self.keep_alive_thread = Thread(group=None, 
                                    target=self.__keep_alive_thread_func, 
                                    name="Thread_Pool keep alive thread")
        self.keep_alive_thread.start()

    def __has_keep_alive_thread(self):

        if self.keep_alive_thread is None:
            return False
        elif self.keep_alive_thread.is_alive() == False:
            return False
            
        return True

    def __keep_alive_thread_func(self):
        self.tasks.join()

    def add_threads_if_needed(self):
        self.remove_finished_threads()
        num_active_threads = len(self.threads)
        num_threads_created = 0
        while num_active_threads < self._max_threads:
            if not self.tasks.empty():
                t = self.add_worker_thread()
                assert(isinstance(t, Thread))
                self.threads.append(t)
                num_active_threads += 1
                num_threads_created += 1
            else:
                break

        if not self.__has_keep_alive_thread():
            self.__start_keep_alive_thread()

    def remove_finished_threads(self):
        while not self.dead_thread_queue.empty():
            try:
                t = self.dead_thread_queue.get_nowait()
                if t is None:
                    break
                else:  
                    for i in range(len(self.threads)-1, 0,-1):
                        if t == self.threads[i]:
                            del self.threads[i]
            except queue.Empty as e:
                return 
        return

    def wait_completion(self):
        """
        Wait for completion of all the tasks in the queue
        """
        self.tasks.join() 

#    @abstractmethod
#    def getNextTask(self):
#        """
#        Iterate over the task queue
#        """
