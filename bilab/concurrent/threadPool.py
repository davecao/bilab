# -*- coding: utf-8 -*-

import time
import traceback

from time import sleep
from threading import Thread, Condition, Event, Lock
try:
    from Queue import Queue, Empty, Full
except ImportError:
    from queue import Queue, Empty, Full

from bilab.concurrent import AbstractLocalThreadPool, ThreadTask, ThreadWorker

__all__ = ['ThreadPool']


class ThreadPool(AbstractLocalThreadPool):
    """
    Pool of threads consuming tasks from a queue
    """
    def __init__(self, num_threads, CheckInterval = 0.5):
        
        super(ThreadPool, self).__init__(num_threads=num_threads,  
                                        WorkerCheckInterval=CheckInterval)

        self._next_thread_id = 0
        #return
        #self.tasks = Queue()
        #self.resultsQueue = Queue()

#        self.__threads = []
#        self.__isJoining = False
#        self.__resizeLock = Condition(Lock())
#        self.__taskLock = Condition(Lock())
#
#        self.setThreadCount(num_threads)

#    def setThreadCount(self, num_threads):
#        """
#        Set the current pool size. Acquires the resizing lock,
#        then calls the internal version to do real work
#        """
#        if self.__isJoining:
#            return False
#
#        self.__resizeLock.acquire()
#
#        try:
#            self.__setThreadCountNolock(num_threads)
#        finally:
#            self.__resizeLock.release()
#
#        return True
#
#    def __setThreadCountNolock(self, num_threads):
#        """
#        Set the current pool size, spawning or terminating threads
#        if necessary. Internal use only; assume the resizing lock is held.
#        """
#
#        while num_threads > len(self.__threads):
#            newThread = Worker(self)
#            self.__threads.append(newThread)
#            #newThread.start()
#
#        # if need to shrink the pool, do so
#        while num_threads < len(self.__threads):
#            self.__threads[0].dismiss()
#            del self.__threads[0]

#    def getThreadCount(self):
#        """ Return the number of threads in the pool """
#        self.__resizeLock.acquire()
#        try:
#            return len(self.__threads)
#        finally:
#            self.__resizeLock.release()

#    def dismissWokers(self, num_threads):
#        """
#        Tell num_threads worker threads to to quit when they're done.
#        """
#        for i in range(min(num_threads, self.__threads)):
#            worker = self.workers.pop()
#            worker.dismiss()
#    def __keep_alive_thread_func(self):
#        self.tasks.join()

    def add_worker_thread(self):
        """
        Implemented by derived class and return a thread object
        """
        w = ThreadWorker(self, self.dead_thread_queue, 
                         self.shutdown_event,
                         self.WorkerCheckInterval)
        w.name = "Thread pool #{}".format(self._next_thread_id)
        self._next_thread_id += 1
        return w

    def add_task(self, taskid, func, *args, **kwargs):
        #def add_task(self, func, *args, **kargs):
        """
        Add a task to the queue
        """
        #callback = kwargs.pop("callback", None)
        #self.tasks.put((func, args, kargs))
#        if self.__isJoining == True:
#            return False
#
#        if not callable(task):
#            return False
#
#        self.__taskLock.acquire()
#        try:
#            #self.tasks.put((taskid, task, args, kwargs, callback))
#            th_task = ThreadTask(task, args, kwargs)
#            self.tasks.put((taskid, th_task))
#            return True
#        finally:
#            self.__taskLock.release()
        th_task = ThreadTask(taskid, func, *args, **kwargs)
        self.tasks.put(th_task)
        self.add_threads_if_needed()

        return th_task

#    def getNextTask(self):
#        """ 
#        Retrieve the next task from the task queue. For use only by Worker 
#        objects contained in the pool. 
#        """
#        self.__taskLock.acquire()
#        try:
#            if self.tasks.empty():
#                #return (None, None, None, None, None)
#                return (None, None)
#            else:
#                return self.tasks.get(0)
#        finally:
#            self.__taskLock.release()

#    def wait_completion(self, waitForTasks=True, waitForThreads = True):
#        """
#        Wait for completion of all the tasks in the queue
#        """
#        # block task queueing
#        self.__isJoining = True
#        # wait for tasks to finish
#        if waitForTasks:
#            while not self.tasks.empty():
#                #print("Current size of tasks:{}".format(self.tasks.qsize()))
#                sleep(.1)
#        #quit for all threads
#        self.__resizeLock.acquire()
#        try:
#            #Wait until all threads have exited
#            if waitForThreads:
#                for th in self.__threads:
#                    th.dismiss()
#                for th in self.__threads:
#                    th.join()
#                    del th
#            self.__setThreadCountNolock(0)
#            # Reset the pool for potential reuse
#            self.__isJoining = False
#        finally:
#            self.__resizeLock.release()
