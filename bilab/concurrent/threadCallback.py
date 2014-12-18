# -*- coding: utf-8 -*-

from threading import Thread, Condition, Event, Lock

__all__ = ['CallBackResults', 'ThreadCallback']


class ThreadCallback(object):
    """
    docstring for ThreadCallback
    """
    lock = Lock()
    def __init__(self):
        super(ThreadCallback, self).__init__()

class CallBackResults(ThreadCallback):
    """
    Store results in callBack function
    of multiple threading function
    """
    def __init__(self):
        super(CallBackResults, self).__init__()
