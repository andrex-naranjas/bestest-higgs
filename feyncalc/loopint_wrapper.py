"""
---------------------------------------------------------------
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          T. Cisneros-Perez
---------------------------------------------------------------
"""
import ctypes
from ctypes import CDLL
import os


class loopint(object):
    """
    Python wrapper for the LoopIntegral C++ shared library
    """
    def __init__(self, workpath="."):
        self.workpath = os.path.dirname(os.path.realpath(__file__))
        
    def loop_integral(self):
        """
        Method to convert the python variables to c++ objects
        """
        # test = ctypes.c_double(MA_val)
        m_lib = ctypes.CDLL(os.path.join(self.workpath+"/LoopIntegral", 'libloopintegral.so'))
        m_lib.loopint_execute.restype = ctypes.c_double
        m_lib.loopint_execute.argtypes = [ctypes.c_double]
        loop_int_value = m_lib.loopint_execute()
        
        return loop_int_value
