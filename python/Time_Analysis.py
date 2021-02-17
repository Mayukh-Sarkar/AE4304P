# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 17:22:09 2021

@author: sarka
"""
import scipy as sc

from stability_analysis import V , A2 , B ,C, D,A_r,B_r,C_r,D_r,b
from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import lsim,lti
from control import matlab as cnt

system_controlled_1 = sc.signal.StateSpace(A2, B, C, D)
system_reduced  = sc.signal.StateSpace(A_r,B_r,C_r,D_r)


dt = 0.05
t = np.linspace(0, 60, num=1201)
u1 = np.zeros((len(t),5))
u2 = np.zeros((len(t),5))
u3 = np.zeros((len(t),5))

N = len(t)
u_g = np.random.randn(1,N)/np.sqrt(dt)  # sqrt(dt) because of lsim characteristics
v_g = np.random.randn(1,N)/np.sqrt(dt);
w_g = np.random.randn(1,N)/np.sqrt(dt);
u1[:,2] = u_g
u2[:,3] = v_g
u3[:,4] = w_g

 


t,y1,x1 = sc.signal.lsim(system_controlled_1, u1, t)
t,y2,x2 = sc.signal.lsim(system_controlled_1, u2, t)
t,y3,x3 = sc.signal.lsim(system_controlled_1, u3, t)
yt = y1+y2+y3
xt = x1 +x2+x3

ay =  V*(xt[:,1] +xt[:,0])


y1r = cnt.lsim(system_reduced, u1, t)
t,y2r,x2r = sc.signal.lsim(A_r,B_r,C_r,D_r, u2, t)
t,y3r,x3r = sc.signal.lsim(system_reduced, u3, t)

plt.plot(t,y1r[:,0])
plt.show()

#xtr = x1r+x2r+x3r

#ayr =  V*(xtr[:,1] +xtr[:,0])



# t = t.reshape(max(t.shape),1)
# ay = ay.reshape(max(ay.shape),1)
# ayr = ayr.reshape(max(ay.shape),1)
# assert t.shape[0] == ay.shape[0]
# assert t.shape[0] == ayr.shape[0]
# plt.plot(t,ay)
# plt.plot(t, yt[:,0])
# #plt.plot(t, y3)
# plt.grid(alpha=0.3)
# plt.xlabel('t')
# plt.show()
# print(ay)

# plt.figure(1)
# plt.plot(t,xt[:,0],'g')
# plt.xlabel('Time ,s')
# plt.ylabel('beta,rad')
# plt.show()

# plt.figure(2)
# plt.plot(t,xtr[:,0],'b')
# plt.xlabel('Time ,s')
# plt.ylabel('beta,rad')
# plt.show()

# plt.figure(3)
# plt.plot(t,xt[:,1],'g')
# plt.xlabel('Time ,s')
# plt.ylabel('phi,rad')
# plt.show()

# plt.figure(4)
# plt.plot(t,xtr[:,1],'b')
# plt.xlabel('Time ,s')
# plt.ylabel('phi,rad')
# plt.show()

# plt.figure(5)
# plt.plot(t,xt[:,2],'g')
# plt.xlabel('Time ,s')
# plt.ylabel('pb/2V,rad')
# plt.show()

# plt.figure(6)
# plt.plot(t,xtr[:,2],'b')
# plt.xlabel('Time ,s')
# plt.ylabel('pb/2V,rad')
# plt.show()

# plt.figure(7)
# plt.plot(t,xt[:,3],'g')
# plt.xlabel('Time ,s')
# plt.ylabel('rb/2V,rad')
# plt.show()

# plt.figure(8)
# plt.plot(t,xtr[:,3],'b')
# plt.xlabel('Time ,s')
# plt.ylabel('rb/2V,rad')
# plt.show()




























