"""
---------------------------------------------------------------
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          T. Cisneros-Perez
---------------------------------------------------------------
"""
import ctypes
from ctypes import (CDLL, POINTER, ARRAY, c_void_p,
                    c_int, byref, c_double, c_char, c_float,
                    c_char_p, create_string_buffer, Structure)
import os
from numpy.ctypeslib import ndpointer


class Complex(Structure):
    """
    Type of object to store the results of the LoopTools functions
    as complex numbers
    """
    _fields_ = [("real", c_double), ("imag", c_double)]
    

class loopint(object):
    """
    Python wrapper for the LoopTools Fortran shared library
    """
    def __init__(self, workpath="."):
        self.workpath = os.path.dirname(os.path.realpath(__file__))
        self.m_lib = ctypes.CDLL(os.path.join(self.workpath+"/LoopTools", 'liblooptools.so'))
        self.m_lib.ltini_()
        
        
    def C0(self, p1s, p2s, p3s, m1, m2, m3):
        """
        Method to obtain C0 values
        """
        m_lib = self.m_lib
        m_lib.c0_.restype = Complex
        valC = m_lib.c0_(byref(c_double(p1s)), byref(c_double(p2s)), byref(c_double(p3s)),
                         byref(c_double(m1)), byref(c_double(m2)), byref(c_double(m3)))
        return complex(valC.real, valC.imag)
