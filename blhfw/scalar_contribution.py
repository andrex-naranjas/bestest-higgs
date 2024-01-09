"""
---------------------------------------------------------------
  Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
           T. Cisneros-Perez
---------------------------------------------------------------
"""
import numpy as np
from blhfw.chromo_moments import ChromoMomentsEM
import scalar_vertex as sv
import masses_blh as mlh

magnetic_electric = ChromoMomentsEM()
qq  = pow(91, 2)
mfi = 173
ms  = 125
mfj = 283.35689658565644

# c_magnetic_scalar = magnetic_electric.chromo_magnetic_scalar(mfi, qq, mfj, ms)
# c_electric_scalar = magnetic_electric.chromo_electric_scalar(mfi, qq, mfj, ms)
# c_magnetic_vector = magnetic_electric.chromo_magnetic_vector(mfi, qq, mfj, ms)
# c_electric_vector = magnetic_electric.chromo_electric_vector(mfi, qq, mfj, ms)
# c_magnetic_photon = magnetic_electric.chromo_magnetic_photon(mfi, qq, mfj, ms)
c_electric_photon = magnetic_electric.chromo_electric_photon(mfi, qq, mfj, ms)


L_scalar = sv.scalar_lagrangian()
vertex_expr = sv.scalar_vertex(L_scalar, "T6bar", "A0", "t")
print("Ready lagrangian")
input()

for _ in range(1000):
    vertex_value = sv.scalar_vertex_value(vertex_expr, beta_value=0, alfa_value=np.pi/2, f_value=1000, y1_value=0.9, y2_value=0.7, y3_value=0.33)
    print(vertex_value, type(vertex_value))

    

f = 1000
beta = 0.55    
mA0, mH0, y3, psi = mlh.EstM41(beta, f)
print("mA0 = ", mA0)
print("mH0 = ", mH0)
print("y3 = ",  y3)
print("psi = ", psi)

mA0 =  140.641824775663
mH0 =  874.426444079677
y3 =  0.968755952873334
psi =  0.169051618495053

MA0 = mA0
MHmm= MA0
MH0 = mH0
Y3  = y3 
Psi = psi

Pi = 3.141592

F = 5000
Lambda = 4*Pi*f
mass_test = mlh.masses_MLH(f, F)
print(mass_test["MA0"])   
print(mass_test["MH0"])   
print(mass_test["Mh0"])   
print(mass_test["mSigma"])
print(mass_test["MHmm"])  
print(mass_test["MPhi"])  
print(mass_test["MPhiMm"])
print(mass_test["MEta0"]) 
print(mass_test["MEtaMm"])
print(mass_test["MZ"])    
print(mass_test["MZp"])   
print(mass_test["MW"])    
print(mass_test["MWp"])   
print(mass_test["MT"])    
print(mass_test["M5"])    
print(mass_test["M6"])    
print(mass_test["M23"])   
print(mass_test["M53"])   
print(mass_test["MB"])
