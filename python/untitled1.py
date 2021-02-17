# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 11:07:28 2021

@author: sarka
"""
from stability_analysis import *
from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import lsim,lti

system_controlled_1 = lti(A2, B, C, D)


dt = 0.05
t = np.linspace(0, 60, num=1201)
u = np.zeros((len(t),5))

print(u)
u_g = np.random.randn(1,len(t))/np.sqrt(dt)
u[:,3] = u_g
#u2 = np.array([nn ,nn ,nn , v_g ,nn])

t,y1,x = lsim(system_controlled_1, u, t)
plt.plot(t, y1[:,0])
plt.grid(alpha=0.3)
plt.xlabel('t')
plt.show()