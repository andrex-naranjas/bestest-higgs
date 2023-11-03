"""
---------------------------------------------------------------
  Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
           T. Cisneros-Perez
---------------------------------------------------------------
"""
from feyncalc.feyncalc_wrapper import loopint
import numpy as np


class LoopIntegral:
    """
    Class that administrates the loop integral calculations done by the C++ class
    The class calls the python wrapper and feeds the functions
    and masses. 
    """
    def __init__(self, bootstrap=False, baryons='', workpath="."):
        self.m_loopint = loopint(workpath)

        
    def loop_integral_value(self, baryons, k_prim, massA, SA_val, L_val, JA_val, SL_val, ModEx_val, bootstrap=False, gamma_val=None, m1=0, m2=0, m3=0):
        """
        Method that calls the wrapper and sums the individual decay widths
        """
        loop_value = self.m_loopint.loop_integral()

