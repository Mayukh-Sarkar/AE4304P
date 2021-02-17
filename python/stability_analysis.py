
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:16:21 2021

@author: sarka
"""
import numpy as np
import scipy as sc
import control as cn
from matplotlib import pyplot as plt
import math

##stability analysis 
V   = 59.9
S   = 24.2
b   = 13.36
mub = 15.5
KX2 = 0.012
KZ2 = 0.037
KXZ = 0.002
CL  = 1.1360

# TURBULENCE PARAMETERS APPROXIMATED POWER SPECTRAL DENSITIES
Lg        = 150 
B         = b/(2*Lg)
sigma     = 2
sigmaug_V = sigma/V
sigmavg   = 1
sigmabg   = sigmavg/V
sigmaag   = sigma/V

Iug0 = 0.0249*sigmaug_V**2
Iag0 = 0.0182*sigmaag**2
tau1 = 0.0991     
tau2 = 0.5545     
tau3 = 0.4159
tau4 = 0.0600     
tau5 = 0.3294     
tau6 = 0.2243

# AIRCRAFT ASYMMETRIC AERODYNAMIC DERIVATIVES 
CYb  =-0.9896     
Clb  =-0.0772     
Cnb  = 0.1638
CYp  =-0.0870     
Clp  =-0.3444     
Cnp  =-0.0108
CYr  = 0.4300     
Clr  = 0.2800     
Cnr  =-0.1930
CYda = 0.0000     
Clda =-0.2349     
Cnda = 0.0286
CYdr = 0.3037     
Cldr = 0.0286     
Cndr =-0.1261
 
                   
Clpw = 0.8*Clp    
Cnpw = 0.9*Cnp
                   
Clrw = 0.7*Clr   
Cnrw = 0.2*Cnr
CYfb = 0
Clfb = 0
Cnfb = 0

#CYfbg = CYfb+0.5*CYr
#Clfbg = Clfb+0.5*Clr
#Cnfbg = Cnfb+0.5*Cnr

# CALCULATION OF AIRCRAFT ASYMMETRIC STABILITY DERIVATIVES
yb   = (V/b)*CYb/(2*mub)
yphi = (V/b)*CL/(2*mub)
yp   = (V/b)*CYp/(2*mub)
yr   = (V/b)*(CYr-4*mub)/(2*mub)
ybg  = yb
ydr  = (V/b)*CYdr/(2*mub)
den  = b*4*mub*(KX2*KZ2-KXZ**2)/V
lb   = (Clb*KZ2+Cnb*KXZ)/den
lp   = (Clp*KZ2+Cnp*KXZ)/den
lr   = (Clr*KZ2+Cnr*KXZ)/den
lda  = (Clda*KZ2+Cnda*KXZ)/den
ldr  = (Cldr*KZ2+Cndr*KXZ)/den
lug  = (-Clrw*KZ2-Cnrw*KXZ)/den
lbg  = lb
lag  = (Clpw*KZ2+Cnpw*KXZ)/den
nb   = (Clb*KXZ+Cnb*KX2)/den
np1   = (Clp*KXZ+Cnp*KX2)/den
nr   = (Clr*KXZ+Cnr*KX2)/den
nda  = (Clda*KXZ+Cnda*KX2)/den
ndr  = (Cldr*KXZ+Cndr*KX2)/den
nug  = (-Clrw*KXZ-Cnrw*KX2)/den
nbg  = nb
nag  = (Clpw*KXZ+Cnpw*KX2)/den
aug1 =-(V/Lg)**2*(1/(tau1*tau2))
aug2 =-(tau1+tau2)*(V/Lg)/(tau1*tau2)
aag1 =-(V/Lg)**2*(1/(tau4*tau5))
aag2 =-(tau4+tau5)*(V/Lg)/(tau4*tau5)
abg1 =-(V/Lg)**2
abg2 =-2*(V/Lg)
bug1 = tau3*math.sqrt(Iug0*V/Lg)/(tau1*tau2)
bug2 = (1-tau3*(tau1+tau2)/(tau1*tau2))*math.sqrt(Iug0*(V/Lg)**3)/(tau1*tau2)
bag1 = tau6*math.sqrt(Iag0*V/Lg)/(tau4*tau5)
bag2 = (1-tau6*(tau4+tau5)/(tau4*tau5))*math.sqrt(Iag0*(V/Lg)**3)/(tau4*tau5)
bbg1 = sigmabg*math.sqrt(3*V/Lg)
bbg2 = (1-2*math.sqrt(3))*sigmabg*math.sqrt((V/Lg)**3)

# STATE- AND INPUT MATRICES
A = np.matrix([[yb ,yphi, yp ,   yr, 0 ,   0,    0 ,   0 ,   ybg , 0],
[0 , 0 ,   2*V/b, 0 , 0   , 0   , 0  ,  0   , 0  ,  0],
[lb, 0 ,   lp   , lr, lug , 0   , lag,  0   , lbg,  0],
[nb, 0 ,   np1   , nr ,nug , 0   , nag,  0   , nbg,  0],
[0 , 0 ,   0    , 0,  0   , 1   , 0  ,  0   , 0  ,  0],
[0 , 0 ,   0    , 0,  aug1, aug2, 0  ,  0   , 0  ,  0],
[0 , 0 ,   0    , 0,  0   , 0 ,   0  ,  1   , 0  ,  0],
[0 , 0 ,   0    , 0,  0   , 0 ,   aag1, aag2, 0  ,  0],
[0 , 0 ,   0    , 0,  0   , 0 ,   0   , 0,    0  ,  1],
[0 , 0 ,   0    , 0,  0   , 0 ,   0   , 0,    abg1 ,abg2]])

B = np.matrix([[0   ,ydr ,0    ,0    ,0],
      [0   ,0   ,0    ,0    ,0],
      [lda ,ldr ,0    ,0    ,0],
      [nda ,ndr ,0    ,0    ,0],
      [0   ,0   ,bug1 ,0    ,0],
      [0   ,0   ,bug2 ,0    ,0],
      [0   ,0   ,0    ,bag1 ,0],
      [0   ,0   ,0    ,bag2 ,0],
      [0   ,0   ,0    ,0    ,bbg1],
      [0   ,0   ,0    ,0    ,bbg2]])

C = np.ones((10,10))
D = np.zeros((10,5))

# SHOW np.linalg.eigENVALUES OF THE UNCONTROLLED SYSTEM
eig_A = np.linalg.eigvals(A)
print('uncontrolled eigen values = ' +  str(eig_A))
sys_uncontrolled = cn.StateSpace(A,B,C,D)
plt.figure(1)
cn.pzmap(sys_uncontrolled,Plot = True,grid =True,title = 'Pole Zero Map of uncontrolled system')
plt.show()


# check for yourself that the spiral mode is not stable, the
# corresponding pole lies in the right-half plane (s = 0.0764).

# THE CESSNA CITATION CE-500 IS NOT STABLE IN SPIRAL MODE (FOR THE cit2a.m 
# FLIGHT CONDITION), HENCE THE FEEDBACK CONTROLLER TO THE AILERON FOR PHI IS 
# USED AS IN : 
#
#   delta_a = K_phi*phi     :  (K_phi for THIS flight condition)
#
# THEREFORE, CONTROLLED AIRCRAFT SYSTEM MATRICES WILL BE USED FOR RESULTS
#
#      A = At 
#
# NO ALTERATIONS MADE FOR cit1a.m !

# NOTE: SPIRAL MODE IS JUST STABLE WITH K_phi ENTERED BELOW
A_r = np.matrix([[0   ,-2*V/b ,0    ,0    ,0    ,0    ,0  ,0],
      [nb    ,nr ,nug  ,0    ,nag  ,0    ,nbg  ,0],
      [0     ,0  ,0    ,1    ,0    ,0    ,0    ,0],
      [0     ,0  ,aug1 ,aug2 ,0    ,0    ,0    ,0],
      [0     ,0  ,0    ,0    ,0    ,1    ,0    ,0],
      [0     ,0  ,0    ,0    ,aag1 ,aag2 ,0    ,0],
      [0     ,0  ,0    ,0    ,0    ,0    ,0    ,1],
      [0     ,0  ,0    ,0    ,0    ,0    ,abg1 ,abg2]])

B_r =np.matrix([[0 , 0   ,0   , 0    ,0],
      [0   ,0  , 0    ,0    ,0],
      [0   ,0  , bug1 ,0   , 0],
      [0   ,0  , bug2 ,0   , 0],
      [0   ,0  , 0    ,bag1, 0],
      [0   ,0  , 0    ,bag2, 0],
      [0   ,0  , 0    ,0   , bbg1],
      [0   ,0  , 0    ,0   , bbg2]])
C_r = np.ones((8,8))
D_r = np.zeros((8,5))

# SHOW np.linalg.eigENVALUES OF THIS CONTROLLED SYSTEM
eig_Ar = np.linalg.eigvals(A_r)
print('Reduced eigen values = ' + str(eig_Ar))
sys_reduced_c = cn.StateSpace(A_r,B_r,C_r,D_r)
plt.figure(2)
cn.pzmap(sys_reduced_c,Plot = True,grid =True,title = 'Pole Zero Map of reduced system')
plt.show()

# FOR EFFECT OF K_phi ON AIRCRAFT RESPONSES, ONE OTHER K_phi IS USED
Kphi = -0.2
K   = np.array([0 ,Kphi, 0 ,0 , 0, 0,  0,  0,  0,  0])
A2 = A-B[:,0]*K

# SHOW np.linalg.eigENVALUES OF THIS CONTROLLED SYSTEM
eig_A_c = np.linalg.eigvals(A2)
print('controlled eig values = ' + str(eig_A_c))
sys_controlled = cn.StateSpace(A2,B,C,D)
plt.figure(3)
cn.pzmap(sys_controlled,Plot = True,grid =True,title = 'Pole Zero Map of controlled system')
plt.show()

print(A)
print(B)
print(A_r)
print(B_r)
print(A2)

