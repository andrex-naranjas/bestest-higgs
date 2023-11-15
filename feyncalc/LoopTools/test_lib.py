import ctypes
import numpy as np
from ctypes import (CDLL, POINTER, ARRAY, c_void_p,
                    c_int, byref,c_double, c_char,c_float,
                    c_char_p, create_string_buffer,Structure)

from numpy.ctypeslib import ndpointer



lib = ctypes.CDLL("./liblooptools.so")



c_int_p = POINTER(c_int)

class Complex(Structure):
    _fields_ = [("real", c_double), ("imag", c_double)]



# def C0(p1s, p2s, p3s, m1, m2, m3):
#     res = Complex()
#     lib.c0_(byref(res),byref(c_double(p1s)),byref(c_double(p2s)),byref(c_double(p3s)),byref(c_double(m1)),byref(c_double(m2)),byref(c_double(m3)))
#     return complex(res.real, res.imag)


lib.ltini_()

p1s = 100
p2s = 2
p3s = 2
m1 = 2
m2 = 2
m3 = 2
res = Complex()

lib.c0_.restype = Complex
lib.b0_.restype = Complex

valC = lib.c0_(byref(c_double(p1s)), byref(c_double(p2s)), byref(c_double(p3s)), byref(c_double(m1)), byref(c_double(m2)), byref(c_double(m3)))
valB = lib.b0_(byref(c_double(p1s)), byref(c_double(p2s)), byref(c_double(p3s)))
print(valB.real, valB.imag, "C0 real e imaginario")
print(valC.real, valC.imag, "B0 real e imaginario")



