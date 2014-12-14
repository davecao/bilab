# -*- coding: utf-8 -*-

from threading import Thread, Condition, Event, Lock

__all__ = ['ThreadWorker']

class ThreadWorker(Thread):
    """
    Thread executing tasks from a given tasks queue
    """
#    threadSleepTime = 0.5

    def __init__(self, tasks, dead_thread_queue, 
                       thread_event, wait_time, **kwargs):

        super(ThreadWorker, self).__init__(**kwargs)

        self.tasks = tasks
        self.dead_thread_queue = dead_thread_queue
        #self._dismissed = Event()
        self.__dismissed = thread_event
        self.wait_time = wait_time

        self.daemon = True
        self.start()
    
    def run(self):
        """
        Repeatedly process the job queue until told to exit.
        """
        #while not self._dismissed.isSet():
        while True:
            try:
                # Thread blocks here, if queue is empty
                #taskid, func, args, kwargs, callback = self.__pool.getNextTask()
                #taskid, th_task = self.__pool.getNextTask()
                taskid, th_task = self.tasks.get(True, self.wait_time)
                print("{}, Task id:{}".format(self.getName(), taskid))
                #try: 
                #    func(*args, **kargs)
                #except Exception, e: 
                #    print e
                #self.tasks.task_done()
                
                # record time
                #task_start_t = time.time()
                #if func is None:
                #    sleep(self.wait_time)
                #elif callback is None:
                #    func(*args, **kwargs)
                #else:
                #    callback(func(*args, **kwargs))
                th_task.run()
            except:
                # Check event
                if self.__dismissed.isSet():
                    # Queue is empty
                    self.dead_thread_queue.put(self)
                    return
                else:
                    # idle state
                    self.dead_thread_queue.put(self)
                    return
            # run the task
            #th_task.run()
            JobsQueued = self.tasks.qsize()
            if JobsQueued > 0:
                JobQText = "Jobs Queued: " + str(self.tasks.qsize())
                JobQText = ('\b' * 40) + JobQText + (' ' * (40 - len(JobQText)))
                print JobQText
            self.tasks.task_done()

    def dismiss(self):
        """
        Set a flag to tell the thread to exit when done with current job.
        """
        self._dismissed.set()
