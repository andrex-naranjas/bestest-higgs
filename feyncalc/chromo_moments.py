"""
---------------------------------------------------------------
  Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
           T. Cisneros-Perez
---------------------------------------------------------------
"""
import numpy as np
from feyncalc.LoopTools.loop_integral import LoopIntegral


class ChromoMomentsEM:
    """
    Class that administrates the loop integral calculations done by the C++ class
    The class calls the python wrapper and feeds the functions
    and masses.
    """
    def __init__(self, workpath="."):
        self.m_loopint = LoopIntegral(workpath)


    def set_A0_values(self, p1s):
        """
        Method to set A0 values
        """
        self.m_A0 = self.m_loopint.A0(p1s)


    def set_B0_values(self, p1s, p2s, p3s):
        """
        Method to set B0 values
        """
        self.m_B0 = self.m_loopint.B0(p1s, p2s, p3s)
        

    def set_C0_values(self, p1s, p2s, p3s, m1, m2, m3):
        """
        Method to set C0 values
        """
        self.m_C0 = self.m_loopint.C0(p1s, p2s, p3s, m1, m2, m3)
        

    def chromo_magnetic(mfi, qq, mfj, ms ):

        Pi = 3.141592
        gs = 1
        B01 = "xx"
        B02 = "xx"
        C01 = "xx"
        A01 = "xx"
        A01 = "xx"

        fm1 = (gs*(4*mfi^4 + 2*B01*mfi^4 - 2*B02*mfi^4 + 2*C01*mfi^6 + A01*(4*mfi^2 - qq) + B02*mfi^2*(6*mfj^2 - 6*ms^2 - qq) + 2*C01*mfi^4*(2*(mfj^2 + ms^2) - qq) - mfi^2*qq + B01*(mfj - ms)*(mfj + ms)*qq + A02*(-4*mfi^2 + qq) +  B01*mfi^2*(-10*mfj^2 + 10*ms^2 + qq) - 2*C01*mfi^2*(3*(mfj^2 - ms^2)^2 - (mfj^2 - 2*ms^2)*qq))*(P1*P2 - S1*S2))




        fm2 = (16*mfi*Pi^2*(-4*mfi^2 + qq)^2) + (gs*mfj*(B01 - B02 - C01*(mfi^2 - mfj^2 + ms^2))*(P1*P2 + S1*S2))/ (8*Pi^2*(4*mfi^2 - qq))


        

        fm = fm1 / fm2
        



        
p1s = 10
p2s = 2
p3s = 2
m1 = 2
m2 = 2
m3 = 2
integral = LoopIntegral()
input()

for _ in range(100000):
    value_integral = integral.A0(p1s)
    print(value_integral,  "   A0")
    value_integral = integral.B0(p1s, p2s, m1)
    print(value_integral,  "   B0")
    value_integral = integral.C0(p1s, p2s, p3s, m1, m2, m3)
    print(value_integral,  "   C0")
