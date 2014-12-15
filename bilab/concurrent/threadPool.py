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

    def __keep_alive_thread_func(self):
        self.tasks.join()

    def  __start_keep_alive_thread(self):
        """ Keep alive threads """
        self.keep_alive_thread = Thread(group=None, 
                                    target=self.__keep_alive_thread_func, 
                                    name="Thread_Pool keep alive thread")
        self.keep_alive_thread.start()

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
        """
        Add a task to the queue
        """
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
