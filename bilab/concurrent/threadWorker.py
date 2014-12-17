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
                th_task = self.tasks.get(True, self.wait_time)
                #print("{}, Task id:{} TaskSize:{} ".format(
                #    self.getName(), th_task.name, self.tasks.qsize()))
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

    def dismiss(self):
        """
        Set a flag to tell the thread to exit when done with current job.
        """
        self.dismissed.set()
        self.dismissed.wait()

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
