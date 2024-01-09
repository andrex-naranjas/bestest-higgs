"""
---------------------------------------------------------------
  Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
           T. Cisneros-Perez
---------------------------------------------------------------
"""
import numpy as np
import math
from feyncalc.LoopTools.loop_integral import LoopIntegral


class ChromoMomentsEM:
    """
    Class to compute
    """
    def __init__(self, workpath="."):
        self.m_loopint = LoopIntegral(workpath)
        
    def show_errors(self):
        """
        Method to obtain a summary of warnings and errors from LoopTools
        """
        self.m_loopint.show_errors()

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
        
    def chromo_magnetic_scalar(self, mfi=10, qq=1, mfj=1, ms=1):
        """
        Method to calculate the magnetic dipole
        """
        Pi = np.pi
        gs = 0.1181 * 2 * pow(Pi, 0.5)
        A01 = self.get_A0_value(pow(mfj, 2))
        A02 = self.get_A0_value(pow(ms, 2))
        B01 = self.get_B0_value(pow(mfi, 2), pow(mfj, 2), pow(ms, 2))
        B02 = self.get_B0_value(qq, pow(mfj, 2), pow(mfj, 2))
        C01 = self.get_C0_value(pow(mfi, 2), pow(mfi, 2), qq, pow(mfj, 2), pow(ms, 2), pow(mfj, 2))

        #print(C01, "C01")
        C01 = -0.0000106797
        P1 = -0.238605
        P2 = 0.238605
        S1 = 0.238605
        S2 = 0.238605

        hmufi = (mfi/gs)*(gs*(4*mfi**4 + 2*B01*mfi**4 - 2*B02*mfi**4 + 2*C01*mfi**6 + A01*(4*mfi**2 - qq) + B02*mfi**2*(6*mfj**2 - 6*ms**2 - qq) + 2*C01*mfi**4*(2*(mfj**2 + ms**2) - qq) - mfi**2*qq + B01*(mfj - ms)*(mfj + ms)*qq + A02*(-4*mfi**2 + qq) + B01*mfi**2*(-10*mfj**2 + 10*ms**2 + qq) - 2*C01*mfi**2*(3*(mfj**2 - ms**2)**2 - (mfj**2 - 2*ms**2)*qq))*(P1*P2 - S1*S2))/(16*mfi*pi**2*(-4*mfi**2 + qq)**2) + (gs*mfj*(B01 - B02 - C01*(mfi**2 - mfj**2 + ms**2))*(P1*P2 + S1*S2))/(8*pi**2*(4*mfi**2 - qq))
        
        return hmufi
        

    def chromo_electric_scalar(self, mfi=10, qq=1, mfj=1, ms=1):
        '''
        Method to calculate the electric chromo dipole
        '''
        Pi = np.pi
        gs = 0.1181 * 2 * pow(Pi, 0.5)
        A01 = self.get_A0_value(pow(mfj, 2))
        A02 = self.get_A0_value(pow(ms, 2))
        B01 = self.get_B0_value(pow(mfi, 2), pow(mfj, 2), pow(ms, 2))
        B02 = self.get_B0_value(qq, pow(mfj, 2), pow(mfj, 2))        
        C01 = self.get_C0_value(pow(mfi, 2), pow(mfi, 2), qq, pow(mfj, 2), pow(ms, 2), pow(mfj, 2))

        I = 1

        #print(C01, "C01")
        C01 = -0.0000106797
        P1 = -0.238605
        P2 = 0.238605
        S1 = 0.238605
        S2 = 0.238605

        hdefi = (mfi/gs)*((-1/8)*gs*mfj*(B01 - B02 - C01*(mfi**2 - mfj**2 + ms**2))*(P2*S1 + P1*S2))/(pi**2*(4*mfi**2 - qq))
        return hdefi


    def chromo_magnetic_vector(self, mfi=10, qq=1, mfj=1, ms=1):
        '''
        Method to calculate the electric chromo dipole
        '''
        Pi = np.pi
        gs = 0.1181 * 2 * pow(Pi, 0.5)
        A01 = self.get_A0_value(pow(mfj, 2))
        A02 = self.get_A0_value(pow(ms, 2))
        B01 = self.get_B0_value(pow(mfi, 2), pow(mfj, 2), pow(ms, 2))
        B02 = self.get_B0_value(qq, pow(mfj, 2), pow(mfj, 2))        
        C01 = self.get_C0_value(pow(mfi, 2), pow(mfi, 2), qq, pow(mfj, 2), pow(ms, 2), pow(mfj, 2))

        I = 1

        #print(C01, "C01")
        C01 = -0.0000106797
        P1 = -0.238605
        P2 = 0.238605
        S1 = 0.238605
        S2 = 0.238605
        hmufi = -(mfi/gs)*(gs*(A01*(4*mfi**2 - qq)*(A1*A2*(mfi**2 + 2*mfi*mfj + mfj**2 + 2*mv**2) + (mfi**2 - 2*mfi*mfj + mfj**2 + 2*mv**2)*V1*V2) - A02*(4*mfi**2 - qq)*(A1*A2*(mfi**2 + 2*mfi*mfj + mfj**2 + 2*mv**2) + (mfi**2 - 2*mfi*mfj + mfj**2 + 2*mv**2)*V1*V2) + mfi**2*(4*mfi**2 - qq)*(A1*A2*(mfi**2 + 2*mfi*mfj + mfj**2 + 2*mv**2) + (mfi**2 - 2*mfi*mfj + mfj**2 + 2*mv**2)*V1*V2) - B02*mfi*(2*A1*A2*mfi*(mfi**4 + 6*mfi**3*mfj + mfi**2*(6*mfj**2 - 11*mv**2) - 2*mfi*(mfj**3 + 5*mfj*mv**2) - 3*(mfj**4 + mfj**2*mv**2 - 2*mv**4)) + A1*A2*(mfi**3 - 3*mfi*mfj**2 - 2*mfj**3 + 10*mfi*mv**2 + 8*mfj*mv**2)*qq + 2*mfi*(mfi**4 - 6*mfi**3*mfj + mfi**2*(6*mfj**2 - 11*mv**2) + 2*mfi*(mfj**3 + 5*mfj*mv**2) - 3*(mfj**4 + mfj**2*mv**2 - 2*mv**4))*V1*V2 + (mfi**3 - 3*mfi*mfj**2 + 2*mfj**3 + 10*mfi*mv**2 - 8*mfj*mv**2)*qq*V1*V2) + 2*C01*mfi*(A1*A2*(mfi**7 - 2*mfi**6*mfj + mfi**4*mfj*(4*mfj**2 - 16*mv**2 - qq) - mfi**5*(5*mfj**2 + 12*mv**2 + qq) - mfj*qq*(mfj**4 - 5*mfj**2*mv**2 + 4*mv**4 + 2*mv**2*qq) + mfi**3*(7*mfj**4 + 17*mv**4 + 8*mv**2*qq + 2*mfj**2*(-6*mv**2 + qq)) + mfi**2*mfj*(-2*mfj**4 + 10*mv**4 + 9*mv**2*qq + mfj**2*(-8*mv**2 + 2*qq)) - mfi*(3*mfj**6 + mfj**4*qq - 3*mfj**2*(3*mv**4 + 2*mv**2*qq) + 2*mv**2*(3*mv**4 + 4*mv**2*qq + qq**2))) + (mfi**7 + 2*mfi**6*mfj - mfi**5*(5*mfj**2 + 12*mv**2 + qq) + mfi**4*mfj*(-4*mfj**2 + 16*mv**2 + qq) + mfi**2*mfj*(2*mfj**4 - 10*mv**4 + mfj**2*(8*mv**2 - 2*qq) - 9*mv**2*qq) + mfj*qq*(mfj**4 - 5*mfj**2*mv**2 + 4*mv**4 + 2*mv**2*qq) + mfi**3*(7*mfj**4 + 17*mv**4 + 8*mv**2*qq + 2*mfj**2*(-6*mv**2 + qq)) - mfi*(3*mfj**6 + mfj**4*qq - 3*mfj**2*(3*mv**4 + 2*mv**2*qq) + 2*mv**2*(3*mv**4 + 4*mv**2*qq + qq**2)))*V1*V2) + B01*(qq*(A1*A2*(mfi**4 + mfj**4 + 6*mfi*mfj*mv**2 + mfj**2*mv**2 - 2*mv**4 + mfi**2*(-2*mfj**2 + 9*mv**2)) + (mfi**4 + mfj**4 - 6*mfi*mfj*mv**2 + mfj**2*mv**2 - 2*mv**4 + mfi**2*(-2*mfj**2 + 9*mv**2))*V1*V2) + 2*mfi**2*(A1*A2*(mfi**4 + 6*mfi**3*mfj + mfi**2*(4*mfj**2 - 9*mv**2) - 6*mfi*mfj*(mfj**2 + mv**2) - 5*(mfj**4 + mfj**2*mv**2 - 2*mv**4)) + (mfi**4 - 6*mfi**3*mfj + mfi**2*(4*mfj**2 - 9*mv**2) + 6*mfi*mfj*(mfj**2 + mv**2) - 5*(mfj**4 + mfj**2*mv**2 - 2*mv**4))*V1*V2))))/(16*mfi*mv**2*pi**2*(-4*mfi**2 + qq)**2)

        return hmufi
        

    def chromo_electric_vector(self, mfi=10, qq=1, mfj=1, ms=1):
        '''
        Method to calculate the electric chromo dipole
        '''
        Pi = np.pi
        gs = 0.1181 * 2 * pow(Pi, 0.5)
        A01 = self.get_A0_value(pow(mfj, 2))
        A02 = self.get_A0_value(pow(ms, 2))
        B01 = self.get_B0_value(pow(mfi, 2), pow(mfj, 2), pow(ms, 2))
        B02 = self.get_B0_value(qq, pow(mfj, 2), pow(mfj, 2))        
        C01 = self.get_C0_value(pow(mfi, 2), pow(mfi, 2), qq, pow(mfj, 2), pow(ms, 2), pow(mfj, 2))

        I = 1

        #print(C01, "C01")
        C01 = -0.0000106797
        P1 = -0.238605
        P2 = 0.238605
        S1 = 0.238605
        S2 = 0.238605

        hdefi = (mfi/gs)*((1/8)*gs*mfj*(B01*(-mfi**2 + mfj**2 - 4*mv**2) + B02*(mfi**2 - mfj**2 + 4*mv**2) + C01*(mfi**4 + mfj**4 - 5*mfj**2*mv**2 + 4*mv**4 - mfi**2*(2*mfj**2 + 3*mv**2) + 2*mv**2*qq))*(A2*V1 - A1*V2))/(mv**2*pi**2*(4*mfi**2 - qq))
        
        return hdefi


    def chromo_magnetic_photon(self, mfi=10, qq=1, mfj=1, ms=1):
        '''
        Method to calculate the magnetic chromo dipole for the photon
        '''
        Pi = np.pi
        gs = 0.1181 * 2 * pow(Pi, 0.5)
        A01 = self.get_A0_value(pow(mfj, 2))
        A02 = self.get_A0_value(pow(ms, 2))
        B01 = self.get_B0_value(pow(mfi, 2), pow(mfj, 2), pow(ms, 2))
        B02 = self.get_B0_value(qq, pow(mfj, 2), pow(mfj, 2))        
        C01 = self.get_C0_value(pow(mfi, 2), pow(mfi, 2), qq, pow(mfj, 2), pow(ms, 2), pow(mfj, 2))

        C00 = 1 # ScalarC0[0, mfi**2, mfi**2, mfj, mfj, 0]

        I = 1

        #print(C01, "C01")
        C01 = -0.0000106797
        P1 = -0.238605
        P2 = 0.238605
        S1 = 0.238605
        S2 = 0.238605
        
        hmufiq0 = (4*mfi**2*(A1*A2*(3*mfi**2 + 4*mfi*mfj + 2*mfj**2) + (3*mfi**2 - 4*mfi*mfj + 2*mfj**2)*V1*V2) + (mfi**2 - mfj**2)*(A1*A2*(7*mfi**2 + 8*mfi*mfj + 5*mfj**2) + (7*mfi**2 - 8*mfi*mfj + 5*mfj**2)*V1*V2)*np.log(mfj**2/(-mfi**2 + mfj**2)) + mfi**2*(A1*A2*(7*mfi**4 + 8*mfi**3*mfj + 6*mfi**2*mfj**2 + 8*mfi*mfj**3 + 3*mfj**4) + (7*mfi**4 - 8*mfi**3*mfj + 6*mfi**2*mfj**2 - 8*mfi*mfj**3 + 3*mfj**4)*V1*V2)*C00)/(64*mfi**4*pi**2)

        return hmufiq0


    def chromo_electric_photon(self, mfi=10, qq=1, mfj=1, ms=1):
        '''
        Method to calculate the electric chromo dipole for the photon
        '''
        Pi = np.pi
        gs = 0.1181 * 2 * pow(Pi, 0.5)
        A01 = self.get_A0_value(pow(mfj, 2))
        A02 = self.get_A0_value(pow(ms, 2))
        B01 = self.get_B0_value(pow(mfi, 2), pow(mfj, 2), pow(ms, 2))
        B02 = self.get_B0_value(qq, pow(mfj, 2), pow(mfj, 2))        
        C01 = self.get_C0_value(pow(mfi, 2), pow(mfi, 2), qq, pow(mfj, 2), pow(ms, 2), pow(mfj, 2))

        C00 = 1 # ScalarC0[0, mfi^2, mfi^2, mfj, mfj, 0]

        I = 1

        #print(C01, "C01")
        C01 = -0.0000106797
        P1 = -0.238605
        P2 = 0.238605
        S1 = 0.238605
        S2 = 0.238605
        
        # preliminary test
        P1 = 1
        P2 = 1
        S1 = 1
        S2 = 1
        ms = 1
        ms = 1

        mfi = 1.0  # Replace with actual value
        mfj = 2.0  # Replace with actual value
        mv = 3.0  # Replace with actual value
        A01 = 4.0  # Replace with actual value
        A02 = 5.0  # Replace with actual value
        A1 = 6.0  # Replace with actual value
        A2 = 7.0  # Replace with actual value
        V1 = 8.0  # Replace with actual value
        V2 = 9.0  # Replace with actual value
        B01 = 10.0  # Replace with actual value
        B02 = 11.0  # Replace with actual value
        C01 = 12.0  # Replace with actual value
        gs = 13.0  # Replace with actual value
        qq = 14.0  # Replace with actual value
        pi = math.pi
        C00 = 1
        # preliminary test
        
        hdefiq0 = ((-1j/8)*mfj*(A2*V1 - A1*V2)*(2*mfi**2 + (mfi**2 - mfj**2)*np.log(mfj**2/(-mfi**2 + mfj**2)) + mfi**2*(mfi**2 + mfj**2)*C00))/(mfi**3*pi**2)

        return hdefiq0


        # # preliminary test
        # P1 = 1
        # P2 = 1
        # S1 = 1
        # S2 = 1
        # ms = 1
        # ms = 1

        # mfi = 1.0  # Replace with actual value
        # mfj = 2.0  # Replace with actual value
        # mv = 3.0  # Replace with actual value
        # A01 = 4.0  # Replace with actual value
        # A02 = 5.0  # Replace with actual value
        # A1 = 6.0  # Replace with actual value
        # A2 = 7.0  # Replace with actual value
        # V1 = 8.0  # Replace with actual value
        # V2 = 9.0  # Replace with actual value
        # B01 = 10.0  # Replace with actual value
        # B02 = 11.0  # Replace with actual value
        # C01 = 12.0  # Replace with actual value
        # gs = 13.0  # Replace with actual value
        # qq = 14.0  # Replace with actual value
        # pi = math.pi
        # # preliminary test
