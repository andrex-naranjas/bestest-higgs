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


c_magnetic_scalar = magnetic_electric.chromo_magnetic_scalar(mfi, qq, mfj, ms, pp1=1, pp2=1, ss1=1, ss2=1)
print(c_magnetic_scalar, "prueba Tzihue")
input()
# c_electric_scalar = magnetic_electric.chromo_electric_scalar(mfi, qq, mfj, ms)
# c_magnetic_vector = magnetic_electric.chromo_magnetic_vector(mfi, qq, mfj, ms)
# c_electric_vector = magnetic_electric.chromo_electric_vector(mfi, qq, mfj, ms)
# c_magnetic_photon = magnetic_electric.chromo_magnetic_photon(mfi, qq, mfj, ms)
# c_electric_photon = magnetic_electric.chromo_electric_photon(mfi, qq, mfj, ms)



tt = np.array(["t", "T", "T5", "T6", "Tb23", "Tb53", "Bp"])
tb = np.array(["tbar", "Tbar", "T5bar", "T6bar", "Tb23bar", "Tb53bar", "Bbar"])
camp = np.array(["A0", "H0", "h0", "sig", "HM", "Hm", "fi0", "fiM", "fim", "eta0", "etaM", "etam"])

L_scalar = sv.scalar_lagrangian()
# vertex_expr = sv.scalar_vertex(L_scalar, "T6bar", "A0", "t")


# tb[[i]]*camp[[j]]*tt[[k]]
for i in range(camp.size):
    for j in range(tt.size):
        vertex_expr = sv.scalar_vertex(L_scalar, tb[j + 1], camp[i], tt[0])
        #print(
        sv.scalar_vertex_value(vertex_expr, beta_value=0.45, alfa_value=-2.161285356607978, f_value="f", y1_value=0.7, y2_value=0.9, y3_value=3.01008, g5_value="g5", pr_value=0.5, pl_value=-0.5)
        # )
        f = 1000+-10
        print(tb[j + 1], camp[i], tt[0])





# L_scalar = sv.scalar_lagrangian()
# vertex_expr = sv.scalar_vertex(L_scalar, "T6bar", "A0", "t")
# print("Ready lagrangian")
# input()

# for _ in range(1000):
#     vertex_value = sv.scalar_vertex_value(vertex_expr, beta_value=0, alfa_value=np.pi/2, f_value=1000, y1_value=0.9, y2_value=0.7, y3_value=0.33)
#     print(vertex_value, type(vertex_value))

    

# f = 1000
# beta = 0.55    
# mA0, mH0, y3, psi = mlh.EstM41(beta, f)
# print("mA0 = ", mA0)
# print("mH0 = ", mH0)
# print("y3 = ",  y3)
# print("psi = ", psi)

# mA0 =  140.641824775663
# mH0 =  874.426444079677
# y3 =  0.968755952873334
# psi =  0.169051618495053

# MA0 = mA0
# MHmm= MA0
# MH0 = mH0
# Y3  = y3 
# Psi = psi

# Pi = 3.141592

# F = 5000
# Lambda = 4*Pi*f
# mass_test = mlh.masses_MLH(f, F)
# print(mass_test["MA0"])   
# print(mass_test["MH0"])   
# print(mass_test["Mh0"])   
# print(mass_test["mSigma"])
# print(mass_test["MHmm"])  
# print(mass_test["MPhi"])  
# print(mass_test["MPhiMm"])
# print(mass_test["MEta0"]) 
# print(mass_test["MEtaMm"])
# print(mass_test["MZ"])    
# print(mass_test["MZp"])   
# print(mass_test["MW"])    
# print(mass_test["MWp"])   
# print(mass_test["MT"])    
# print(mass_test["M5"])    
# print(mass_test["M6"])    
# print(mass_test["M23"])   
# print(mass_test["M53"])   
# print(mass_test["MB"])
