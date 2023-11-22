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
        

    def chromo_magnetic(self, mfi=10, qq=1, mfj=1, ms=1):
        """
        Method to set C0 value
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

        ff=(gs*(4*pow(mfi, 4) + 2*B01*pow(mfi, 4) - 2*B02*pow(mfi, 4) + 2*C01*pow(mfi, 6) + A01*(4*pow(mfi, 2) - qq) + B02*pow(mfi, 2)*(6*pow(mfj, 2) - 6*pow(ms, 2) - qq) + \
                2*C01*pow(mfi, 4)*(2*(pow(mfj, 2) + pow(ms, 2)) - qq) - pow(mfi, 2)*qq + B01*(mfj - ms)*(mfj + ms)*qq + A02*(-4*pow(mfi, 2) + qq) + \
                B01*pow(mfi, 2)*(-10*pow(mfj, 2) + 10*pow(ms, 2) + qq) - 2*C01*pow(mfi, 2)*(3*pow((pow(mfj, 2) - pow(ms, 2)), 2) - (pow(mfj, 2) - 2*pow(ms, 2))*qq))*(P1*P2 - S1*S2))/ \
        (16*mfi*pow(Pi, 2)*pow((-4*pow(mfi, 2) + qq), 2)) 

        fff = (gs*mfj*(B01 - B02 - C01*(pow(mfi, 2) - pow(mfj, 2) + pow(ms, 2)))*(P1*P2 + S1*S2))/ (8*pow(Pi, 2)*(4*pow(mfi, 2) - qq))
        fm = ff + fff
        
        return (fm * mfi) / gs


def sample_gauss(mu, sigma):
    return np.random.normal(mu, sigma, 10000)

def random(sample):
    return np.random.choice(sample, size=None)
        
magnetic = ChromoMomentsEM()

import time
start = time.time()


qq = pow(91, 2)
mfi = 173
ms = 125
mfj = 283.35689658565644


gauss_qq  = sample_gauss(pow(91, 2), 0.002)
gauss_mfi = sample_gauss(173, 0.01)
gauss_ms  = sample_gauss(125, 0.11)
gauss_mfj = sample_gauss(283.35689658565644, 0.01)
sampled_cm = ([])


for _ in range(10000):
    c_magnetic = magnetic.chromo_magnetic(mfi=random(gauss_mfi),
                                          qq=random(gauss_qq),
                                          mfj=random(gauss_mfj),
                                          ms=random(gauss_ms))
    
    #print(c_magnetic.real)
    sampled_cm = np.append(sampled_cm, c_magnetic.real)

end = time.time()
print(end - start)
magnetic.show_errors() # check for errors and warnings
print(sampled_cm.mean(), sampled_cm.std())

import matplotlib.pyplot as plt


# Calculate mean and standard deviation
mean = np.mean(sampled_cm)
std_dev = np.std(sampled_cm)

# Plotting the Gaussian distribution
plt.hist(sampled_cm, bins=100, density=True, alpha=0.7, color='skyblue')

# Adding legend with mean and standard deviation information
plt.legend([f'Mean: {mean:.8f}\nError: {std_dev:.10f}'], loc='lower left')

# Define ranges for x and y axes
#plt.xlim([0.00003, 0.1])  # Set x-axis range from -3 to 3
#plt.xlim([-100, 100])  # Set x-axis range from -3 to 3
#plt.ylim([0, 0.5])  # Set y-axis range from 0 to 0.5

# Display the plot
plt.xlabel('magnetic moment [a.u.]')
plt.ylabel('PDF')
plt.title('Chromomagnetic dipole')
# Save figure as a PDF file
plt.savefig('chromo_magnetic.pdf', format='pdf')
plt.show()

