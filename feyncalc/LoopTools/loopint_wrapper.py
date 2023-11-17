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
from feyncalc.LoopTools import tlt_coeficients


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
    def __init__(self, workpath=".", mu=1, delta=-1.):
        self.workpath = os.path.dirname(os.path.realpath(__file__))
        self.m_lib = ctypes.CDLL(os.path.join(self.workpath, "liblooptools.so"))
        self.m_lib.ltini_()
        # set regularization parameters
        self.m_lib.setmudim_(byref(c_double(pow(mu ,2))) ) # mu^2
        self.m_lib.setdelta_(byref(c_double( delta )))     # delta
        # self.m_lib.setlambda_(byref(c_double( 0 )))      # lambda
        # self.m_lib.setminmass_(byref(c_double( 0 )))     # min mass
        self.m_lib.getmudim_.restype = c_double
        self.m_lib.getdelta_.restype = c_double
        print("\mu is %lf"  %self.m_lib.getmudim_())
        print("\Delta is %lf " %self.m_lib.getdelta_())

        
    def A0(self, x):
        """
        Method to obtain A0 values
        """
        m_lib = self.m_lib
        m_lib.a0_.restype = Complex
        valA = m_lib.a0_(byref(c_double(x)) )  #, c_int_p(c_int(0)))
        return complex(valA.real, valA.imag)


    def B0(self, ps, m1, m2):
        """
        Method to obtain B0 values
        """
        m_lib = self.m_lib
        m_lib.b0_.restype = Complex
        valB = m_lib.b0_(byref(c_double(ps)), byref(c_double(m1)), byref(c_double(m2)))
        return complex(valB.real, valB.imag)
        
        
    def C0(self, p1s, p2s, p3s, m1, m2, m3):
        """
        Method to obtain C0 values
        """
        m_lib = self.m_lib
        m_lib.c0_.restype = Complex
        valC = m_lib.c0_(byref(c_double(p1s)), byref(c_double(p2s)), byref(c_double(p3s)),
                         byref(c_double(m1)), byref(c_double(m2)), byref(c_double(m3)))
        return complex(valC.real, valC.imag)


