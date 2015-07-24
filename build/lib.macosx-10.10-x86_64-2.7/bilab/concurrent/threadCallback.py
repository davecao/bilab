# -*- coding: utf-8 -*-

from threading import Thread, Condition, Event, Lock, BoundedSemaphore

__all__ = ['CallBackResults', 'ThreadCallback']


class ThreadCallback(object):
    """
    docstring for ThreadCallback
    """
    #lock = Lock()
    def __init__(self, max_conn):
        super(ThreadCallback, self).__init__()
        self.lock = BoundedSemaphore(value=max_conn)

class CallBackResults(ThreadCallback):
    """
    Store results in callBack function
    of multiple threading function
    """
    def __init__(self, max_conn):
        super(CallBackResults, self).__init__(max_conn)
