"""
---------------------------------------------------------------
  Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
           T. Cisneros-Perez
---------------------------------------------------------------
"""
from feyncalc.LoopTools.loopint_wrapper import loopint
import numpy as np


class LoopIntegral:
    """
    Class that administrates the loop integral calculations done by the C++ class
    The class calls the python wrapper and feeds the functions
    and masses. 
    """
    def __init__(self, workpath="."):
        self.m_loopint = loopint(workpath)

        
    def C0(self, p1s, p2s, p3s, m1, m2, m3):
        """
        Method that calls the wrapper and sums the individual decay widths
        """
        loop_value = self.m_loopint.C0(p1s, p2s, p3s, m1, m2, m3)
        return loop_value


p1s = 10
p2s = 2
p3s = 2
m1 = 2
m2 = 2
m3 = 2
integral = LoopIntegral()

for _ in range(100000):
    value_integral = integral.C0(p1s, p2s, p3s, m1, m2, m3)
    print(value_integral)
