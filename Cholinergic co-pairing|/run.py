#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 11:31:06 2021

@author: ines
"""
import numpy as np
import matplotlib.pyplot as plt

import E
import fast_neurons
import OLM

from ode import body_new, noise_generation, periodic_rect

# Define time step and total time of simulation
tmin = 0
tmax = 2400e3
dt = 0.02
steps = (tmax - tmin) / dt  # number of steps
T = np.linspace(tmin, tmax, int(steps))  # time vector
T_min = T * 1e-3 / 60 # time in minutes (for plotting)

# Parameters for neurotransmitters
ACh_par ={'init' : 900, 'dur' : 5, 'per' : 60e3}
Glu_par = {'init': 1000, 'dur': 5,'per': 60e3}

# Fast interneuron initial values
VI = np.zeros(len(T))
V_I=-64
m_I = 0.0
h_I = 0.0
n_I = 0.0
r_A = 0
r_G = 0
GABA_I = 0

# OLM initial values
VO = np.zeros(len(T))
V_O = -60
Ca_a7 = 0
Ca_is = 0.44e-3
h_O = 0.0
n_O = 0.0
m_O = 0.0
p = 0
hf = 0
hs = 0
ra7 = 0
ra_O = 0
GABA_O = 0
GABAO = np.zeros(len(T))

# E-cell initial values
ra_E = 0 
rn_E = 0  
rg_E = 0  
Ca=0 
V_Ed = -67 
g = E.AMPA['g'] #AMPAR maximal conductance


CaA7 = np.zeros(len(T))
CaIS = np.zeros(len(T))

Ca_Ed = np.zeros(len(T))
gA = np.zeros(len(T))

# Loop for file generation
x0 = (V_O, m_O, h_O, n_O, p, hf, hs, ra7, Ca_a7, Ca_is, GABA_O, 
      V_I, r_G, r_A, m_I, h_I, n_I, 
      V_Ed, ra_E, rn_E, rg_E, Ca, g, 0., 0.)


noise_O, noise_I, noise_Ed = noise_generation(T, dt) # Different for each step
ACh = periodic_rect(T,ACh_par['init'],ACh_par['dur'], ACh_par['per'])
Glu = periodic_rect(T,Glu_par['init'],Glu_par['dur'], Glu_par['per'])



file = open('EPSC_v1.txt','w')

for i,time in enumerate(T):
    
    # Define co-pairing period  
    if (530e3 + ACh_par['init'])>time or time>(970e3 + ACh_par['init'] ) : ACh[i] = 0 #tdis = 8mints
   # if (530e3 + ACh_par['init'])>time or time>(790e3+ ACh_par['init']) : ACh[i] = 0 #tdis = 5mints
   
    
    # group arguments to be sent to integrate function
    arg = dt, ACh[i], Glu[i], noise_O[i], noise_I[i], noise_Ed[i]
    
    # integrate ODE's
    dynamics = body_new(x0, arg)
    VO[i] = dynamics[0]
    VI[i] = dynamics[11]
    GABAO[i] = dynamics[10]
    CaA7[i] = dynamics[8]
    CaIS[i] = dynamics[9]
    Ca_Ed[i] = dynamics[-4]
    gA[i] = dynamics[-3]
    
    x0 = dynamics
 
    # save local maximum of EPSC
    if Glu[i-1]>Glu[i]:
        in_E = dynamics[-1]
        ia_E = dynamics[-2]
        
        file.write("%s %s\n" % (T_min[i-1], -ia_E - in_E))
        print(T_min[i-1])
     
file.close()



