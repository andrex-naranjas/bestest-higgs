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

    def show_errors(self):
        self.m_loopint.show_errors()

        
    def A0(self, p1s):
        """
        Method that calls the wrapper LoopTools A0
        """
        loop_value = self.m_loopint.A0(p1s)
        return loop_value


    def B0(self, p1s, p2s, p3s):
        """
        Method that calls the wrapper LoopTools B0
        """
        loop_value = self.m_loopint.B0(p1s, p2s, p3s)
        return loop_value
    

    def C0(self, p1s, p2s, p3s, m1, m2, m3):
        """
        Method that calls the wrapper  LoopTools C0
        """
        loop_value = self.m_loopint.C0(p1s, p2s, p3s, m1, m2, m3)
        return loop_value
