from sympy import solve, sqrt, sin, cos, log, Symbol
#import numpy as np

# first method
Pi = 3.14159265
y1 = 0.7
v = 246
y2 = 0.9
lambda0 = 4*Pi
m4 = 500
gA = 0.9231974770385101
gB = 0.9231974770385101
g = 0.6528
gp = 0.3528
Mh0 = 125

# second method
Mt=172

m4=500
MZ=91
MW=81

def EstM41(beta, f):
    '''
    Method to compute the mass of
    Arguments: beta angle and f symmetry breaking scale
    Returns: mA0, mH0, y3, psi
    '''
    x = Symbol('x') # using x for MA0
    MasaA0 = solve(sqrt((x**2 + v**2 * lambda0)/2 - sqrt((x**2 + v**2 * lambda0)**2/4
                   + v**4 * lambda0**2 * sin(2 * beta)**2 -v**2 * lambda0 *
                   (x**2 + v**2 * lambda0) * sin(2 * beta)**2)) - 125, x)
    mA0 = MasaA0[1]
    mH0 = sqrt((mA0**2 + v**2 * lambda0)/2 + sqrt((mA0**2 + v**2 * lambda0)**2/4
                   + v**4 * lambda0**2 * sin(2 * beta)**2 - v**2 * lambda0 *
                   (mA0**2 + v**2 * lambda0) * sin(2 * beta)**2))

    yt  = 172.76 / (v * sin(beta))
    y   = Symbol('y') # using y for y3
    y33 = solve(((3*y1*y2*y)/(sqrt(y1**2 + y2**2)*sqrt(y1**2 + y**2)))-yt, y)
    y3  = y33[0]
    a   = (27*(f**2))/(8*(Pi**2)*(lambda0)*(cos(beta)**2)*(v**2))
    b   = ((y1**2 * y2**2 * y3**2)/(y2**2 - y3**2)) * log((y1**2 + y2**2)/(y1**2 + y3**2))
    psi = a*b
    return mA0, mH0, y3, psi


def masses_MLH(f, F):
    '''
    Method to compute the mass of
    Arguments: beta angle and f symmetry breaking scale
    Returns: a dictionary
    '''
    MWp    = sqrt( (1/4)*( (gA**2) + (gB**2))*((f**2) + (F**2)) - 81)
    a      = (16*9*(F**2)*(gA**2)*(gB**2))/(3*512*(Pi**2))
    MPhi   = sqrt((( (f**4) + (F**4))*(m4**2))/((F**2)*((f**2) + (F**2))) + a*(log( (Lambda**2)/MWp )) )
    MPhiMm = sqrt((( (f**4) + (f**2)*(F**2) + (F**4))*(m4**2))/((F**2)*((f**2) + (F**2))) + 
                   (3*(F**2)*(gA**2)*(gB**2)*log( (Lambda**2)/MWp))/(32*(Pi**2)))
    MEta0  = m4
    MEtaMm = sqrt((m4**2) + (2*(f**2)*3*(gp**2)*(Lambda**2))/(128*(Pi**2)*(F**2)))
    mSigma = sqrt(2*lambda0) * f
    MZp    = sqrt(0.25*((gA**2) + (gB**2))*((f**2) + (F**2)) - 0.25*(g**2)*(v**2))
    M6     = y1 * f
    M23    = y1 * f
    M53    = y1 * f
    M5     = sqrt(-((Mt**2*((y1**2) + (y2**2)))/((y2**2) - (y3**2))) + (f**2)*((y1**2) + (y3**2)))
    MB     = 1.140175425099138 * f
    MT     = sqrt((f**2)*((y1**2) + (y2**2)) + ((Mt**2)*((y1**2) + (y3**2)))/( (y2**2) - (y3**2)))
    masses = {
        "MA0":    MA0,
        "MH0":    MH0,
        "Mh0":    Mh0,
        "mSigma": mSigma,
        "MHmm":   MHmm,
        "MPhi":   MPhi,
        "MPhiMm": MPhiMm,
        "MEta0":  MEta0,
        "MEtaMm": MEtaMm,
        "MZ":     MZ,
        "MZp":    MZp,
        "MW":     MW,
        "MWp":    MWp,
        "MT":     MT,
        "M5":     M5,
        "M6":     M6,
        "M23":    M23,
        "M53":    M53,
        "MB":     MB}
    return masses
