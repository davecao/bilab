# -*- coding: utf-8 -*-

__all__ = []


# decorator borrowed from Mozilla mxr
def abstractmethod(method):
    if PY3K:
        line = method.__code__.co_firstlineno
        filename = method.__code__.co_filename
    else:
        line = method.func_code.co_firstlineno
        filename = method.func_code.co_filename

    @wraps(method)
    def not_implemented(*args, **kwargs):
        raise NotImplementedError(
            'Abstract method %s at File "%s", line %s'
            'should be implemented by a concrete class' %
            (repr(method), filename, line))
    return not_implemented


class ClassRegistry(type):
    """ Register all subclasses """
    def __init__(cls, name, bases, nmspc):
        super(ClassRegistry, cls).__init__(name, bases, nmspc)
        if not hasattr(cls, 'registry'):
            cls.registry = set()
        cls.registry.add(cls)
        cls.registry -= set(bases)  # Remove base classes

    # Meta methods, called on class objects:
    def __iter__(cls):
        return iter(cls.registry)

    def __str__(cls):
        if cls in cls.registry:
            return cls.__name__
        return cls.__name__ + ":" + ', '.join([sc.__name__ for sc in cls])


class MDFileIOParser(object):
    """
        Base class
    """
    def __init__(self):
        super(MDFileIOParser, self).__init__()

    @abstractmethod
    def parse(self):
        """
        Args:
            s   --  a string

        Returns:
            An object of Universe -- bilab.structure.MDSim.Universe
        """
        pass 

