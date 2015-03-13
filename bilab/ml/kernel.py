# -*- coding: utf-8 -*-

__all__ = ["Kernel"]

class ClassRegistry(type):
    """ Register all subclasses """
    def __init__(cls, name, bases, nmspc):
        super(ClassRegistry, cls).__init__(name, bases, nmspc)
        if not hasattr(cls, 'registry'):
            cls.registry = set()
        cls.registry.add(cls)
        cls.registry -= set(bases) #Remove base classes
    # Meta methods, called on class objects:
    def __iter__(cls):
        return iter(cls.registry)

    def __str__(cls):
        if cls in cls.registry:
            return cls.__name__
        return cls.__name__ + ":" + ', '.join([sc.__name__ for sc in cls])


class Kernel(object):
    """ base class

    .. note::
        This class uses two patterns, composite and registry
    """
    __metaclass__ = ClassRegistry

    def __init__(self, kernel=None, *args, **kwargs):
        """Initialization.

        Args:
           kernel (class):

        Kwargs:
        """
        super(Kernel, self).__init__()
        self.__kernel = None

        if kernel:
            # instance, developped by a user or subclasses of Kernel
            self.__kernel = kernel()
        elif self.__isstr(kernel):
            # is string
            self.__kernel = self.__create(kernel, *args, **kwargs)

    def __create(self, clsname, *args, **kwargs):
        """ Create an instance for a given class in str or name

        """
        obj = None
        for cls in self.registry:
            if clsname == cls.__name__:
                obj = cls(args)
            elif clsname == cls:
                obj = cls(args)
        if obj:
            # Initialize object
            obj.__init__(*args, **kwargs)
        else:
            print("Unknown class {}".format())
        return obj

    def __str__(self):
        #OK:    return self.__kernel.__str__()
        #Failed return self.__getattr__(self.__kernel, '__str__')
        kernel_func_str_ = self.__getattr__(self.__kernel, '__str__')
        return kernel_func_str_()

    def __repr__(self):
        #return self.__kernel.__repr__()
        kernel_func_str_ = self.__getattr__(self.__kernel, '__repr__')
        return kernel_func_str_()

    def __call__(self, *args, **kwargs):
        return self.__kernel(*args, **kwargs)

    def __getattr__(self, obj, attr):
        """ get delegation to the object """
        #return getattr(obj, attr)
        try:
            return self.__kernel.__getattribute__(attr)
        except AttributeError:
            raise AttributeError('{0} object has no attribute `{1}`'
                .format(self.__class__.__name__, attr))

    def __isstr(self, s):
        try:
            return isinstance(s, basestring)
        except NameError:
            return isinstance(s, str)
