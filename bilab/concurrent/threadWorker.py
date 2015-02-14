# -*- coding: utf-8 -*-
import logging
from threading import Thread, Condition, Event, Lock

__all__ = ['ThreadWorker']

logging.basicConfig(level=logging.DEBUG, 
    format='[%(levelname)s] (%(threadName)-10s) %(funcName)s %(message)s',)

class ThreadWorker(Thread):
    """
    Thread executing tasks from a given tasks queue
    """
#    threadSleepTime = 0.5

    def __init__(self, tasks, dead_thread_queue, 
                       thread_event, wait_time, 
                       verbose=False, 
                       **kwargs):

        super(ThreadWorker, self).__init__(**kwargs)

        self.verbose = verbose
        self.tasks = tasks
        self.dead_thread_queue = dead_thread_queue
        self.dismissed = thread_event
        self.wait_time = wait_time
        if wait_time < 0:
            self.getFuncName = {'get':self.get_nowait}
        else:
            self.getFuncName = {'get':self.get}
        self.daemon = True
        self.start()
    
    def run(self):
        """
        Repeatedly process the job queue until told to exit.
        """
        #while not self.dismissed.isSet():
        while True:
            try:
                # Thread blocks here, if queue is empty
                #taskid, func, args, kwargs, callback = self.__pool.getNextTask()
                #taskid, th_task = self.__pool.getNextTask()
                #th_task = self.tasks.get(True, self.wait_time)
                th_task = self.tasks.get_nowait()
                #th_task = self.__getattr__(self.getFuncName['get'])()
#                print("{}, Task id:{} TaskSize:{} ".format(
#                    self.getName(), th_task.name, self.tasks.qsize()))
                #print("Task: {}".format(th_task))
                # run the task
                #if self.verbose:
                #    logging.debug('Starting {}'.format(th_task))
                th_task.run()

                #if self.verbose:
                #    logging.debug('Exiting {}'.format(th_task))
                #    JobsQueued = self.tasks.qsize()
                #    if JobsQueued > 0:
                #        JobQText = "Jobs Queued: " + str(self.tasks.qsize())
                #        #JobQText = ('\b' * 40) + JobQText + (' ' * (40 - len(JobQText)))
                #        logging.debug(JobQText)

                # Unblocks the queue.
                self.tasks.task_done()

            except Exception as e:
                # Check event
                if self.dismissed.isSet():
                    # Queue is empty
                    self.dead_thread_queue.put(self)
                    return
                else:
                    # idle state
                    self.dead_thread_queue.put(self)
                    return
    def get_nowait():
        """
            get a task from task queue immediately
        """
        logging.debug(' get_nowait is called.')
        return self.task.get_nowait()

    def get():
        """
            get a task while blocking for timeout
        """
        logging.debug(' get is called.')
        return self.task.get(True, self.wait_time)

    def dismiss(self):
        """
        Set a flag to tell the thread to exit when done with current job.
        """
        self.dismissed.set()
        self.dismissed.wait()

    def __getattr__(self, attr_name):
        try:
            return self.__implementation.__getattribute__(attr_name)
        except AttributeError:
            raise AttributeError('{0} object has no attribute `{1}`'
                .format(self.__class__.__name__, attr_name))
#    def exit(self):
#        """ 
#        Force to exit 
#        Raise the SystemExit exception
#        """
#        try:
#            self.exit()
#        except SystemExit as e:
#            print("Error: Thread {} is failed to exit".format(self.getName()))
#            error_message = "{}\n".format(traceback.format_exc())
