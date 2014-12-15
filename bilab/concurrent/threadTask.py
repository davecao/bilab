# -*- coding: utf-8 -*-

from functools import  wraps
from threading import Event, Thread, currentThread
import traceback
import math

__all__ = ['Task', 'TaskWithEvent', 'ThreadTask']

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

def time_format(elapsed_time):
    # generate the time elapsed string for output
    seconds = math.fmod(t_delta, 60)
    seconds_str = "%02.5g" % seconds
    time_str = str(time.strftime('%H:%M:', time.gmtime(t_delta))) + seconds_str
    return time_str

class Task(object):
    """
    A base class for a task assigned to a pool.
    """

    def __init__(self, name, func, *args, **kwargs):
        """
        Args:
            name(str): task id
            func(object): the command 
        """
        self.name = name 
        self.func = func
        self.callback = kwargs.pop("callback", None) 
        self.args = args
        self.kwargs = kwargs

    @abstractmethod
    def wait(self):
        """
        Wait for task to complete, does not return a value
        """
    @abstractmethod
    def wait_return(self):
        """
        The task has return values
        """

    @property
    def iscompleted(self):
        """ 
        Non-blocking test to determine if task has completed.  
        No exception is raised if the task raised an exception 
        during execution until wait or wait_return is called.
        
        :return: True if the task is completed, otherwise False
        """

class TaskWithEvent(Task):
    """
    Task object uses threading.Event()
    """
    def __init__(self, name, func, *args, **kwargs):
        super(TaskWithEvent, self).__init__(name, func, *args, **kwargs)
        self.completed = Event()

    @property
    def iscompleted(self):
        return self.completed.isSet()
    
    def wait(self):
        self.completed.wait()

class ThreadTask(TaskWithEvent):
    """
    Class for storing info of a threaded task
    """
    @property
    def exception(self):
        return self._exception

    @exception.setter
    def exception(self, value):
        self._exception = value

    def __init__(self, name, func, *args, **kwargs):

        super(ThreadTask, self).__init__(name, func, *args, **kwargs)
        #self.callback = kwargs.pop("callback", None)
        self.verbose = kwargs.pop("verbose", False)
        self.returned_value = None
        self._exception = None

    def __str__(self):
        """ Representative string """
        msg = self.func.__name__ 
        for i, val in enumerate(self.args):
            msg += 'arg[{}]={} '.format(i, val)
        for name, val in self.kwargs.items():
            msg += 'kwargs[{}]={} '.format(name, val)

        return msg

    def __repr__(self):
        return self.__str__()

    def wait(self):
        #self._iscompleted.wait()
        super(ThreadTask, self).wait()
        if not self.exception is None:
            raise self.exception

    def wait_return(self):
        """
        Waits until the function has completed execution and 
        returns the value returned by the function pointer
        """
        self.wait()
        if not self.exception is None:
            raise self.exception
        return self.returned_value

    def run(self):
        """ Run func """
        # record time
        if self.verbose:
            task_start_t = time.time()
        #if func is None:
        #    sleep(self.wait_time)
        #elif callback is None:
        #    func(*args, **kwargs)
        #else:
        #    callback(func(*args, **kwargs))
        currTh = currentThread()
        #print("Thread -{}-, Task name: -{}-".format(currTh.getName(),
        #                                        self.name))
        try:
            if self.func is None:
                self.returned_value = self.func()
            elif self.callback is None:
                self.returned_value = self.func(*self.args, **self.kwargs)
            else:
                #self.returned_value = self.func(*self.args, **self.kwargs)
                #self.callback(self.returned_value)
                self.callback(self.func(*self.args, **self.kwargs))

        except Exception as e:
            self.exception = e
            error_message = "\n*** {0}\n{1}\n".format(
                self.name, traceback.format_exc())
            print error_message

        if self.verbose:
            task_end_time = time.time()
            t_delta = task_end_time - task_start_time
            print("elapsed time: {}".format(time_format(t_delta)))
        # mark the object event as completed
        self.completed.set()
