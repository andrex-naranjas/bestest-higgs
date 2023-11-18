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


    def get_A0_value(self, p1s):
        """
        Method to set A0 value
        """
        return self.m_loopint.A0(p1s)


    def get_B0_value(self, p1s, p2s, p3s):
        """
        Method to set B0 value
        """
        return self.m_loopint.B0(p1s, p2s, p3s)
        

    def get_C0_value(self, p1s, p2s, p3s, m1, m2, m3):
        """
        Method to set C0 value
        """
        return self.m_loopint.C0(p1s, p2s, p3s, m1, m2, m3)
        

    def chromo_magnetic(self, mfi=10, qq=1, mfj=1, ms=1):

        Pi = 3.141592
        gs = 1
        B01 = self.get_B0_value(pow(mfi, 2), pow(mfj, 2), pow(ms, 2))
        B02 = self.get_B0_value(qq, pow(mfj, 2), pow(mfj, 2))
        C01 = self.get_C0_value(pow(mfi, 2), pow(mfi, 2), qq, pow(mfj, 2), pow(ms, 2), pow(mfj, 2))
        A01 = self.get_A0_value(pow(mfj, 2))
        A02 = self.get_A0_value(pow(ms, 2))

        P1 = -0.238605
        P2 = 0.238605
        S1 = 0.238605
        S2 = 0.238605

        ff=(gs*(4*pow(mfi, 4) + 2*B01*pow(mfi, 4) - 2*B02*pow(mfi, 4) + 2*C01*pow(mfi, 6) + A01*(4*pow(mfi, 2) - qq) + B02*pow(mfi, 2)*(6*pow(mfj, 2) - 6*pow(ms, 2) - qq) + \
                2*C01*pow(mfi, 4)*(2*(pow(mfj, 2) + pow(ms, 2)) - qq) - pow(mfi, 2)*qq + B01*(mfj - ms)*(mfj + ms)*qq + A02*(-4*pow(mfi, 2) + qq) + \
                B01*pow(mfi, 2)*(-10*pow(mfj, 2) + 10*pow(ms, 2) + qq) - 2*C01*pow(mfi, 2)*(3*pow((pow(mfj, 2) - pow(ms, 2)), 2) - (pow(mfj, 2) - 2*pow(ms, 2))*qq))*(P1*P2 - S1*S2))/ \
        (16*mfi*pow(Pi, 2)*pow((-4*pow(mfi, 2) + qq), 2)) 

        fff = (gs*mfj*(B01 - B02 - C01*(pow(mfi, 2) - pow(mfj, 2) + pow(ms, 2)))*(P1*P2 + S1*S2))/ (8*pow(Pi, 2)*(4*pow(mfi, 2) - qq))
        fm = ff + fff
        
        return fm

        
magnetic = ChromoMomentsEM()
input()

import time
start = time.time()


qq = pow(91, 2)
mfi = 173
ms = 125
mfj = 283.35

for _ in range(100):
    c_magnetic = magnetic.chromo_magnetic(mfi=mfi, qq=qq, mfj=mfj, ms=ms)
    print(c_magnetic)

end = time.time()
print(end - start)
    # value_integral = integral.A0(p1s)
    # print(value_integral,  "   A0")
    # value_integral = integral.B0(p1s, p2s, m1)
    # print(value_integral,  "   B0")
    # value_integral = integral.C0(p1s, p2s, p3s, m1, m2, m3)
    # print(value_integral,  "   C0")
