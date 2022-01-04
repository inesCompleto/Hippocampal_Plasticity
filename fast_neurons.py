"""
Functions associated with dynamics of fast-spiking interneuron and 
                                      pyr neuron
"""

import numpy as np
        
# fast-spiking interneuron
def dhdt(h, v):
    a_h = 0.128*np.exp(-(v+50)/18)
    b_h = 4/(np.exp(-(v+27)/5)+1)
    h_inf = a_h/(a_h + b_h)
    tau_h = 1/(a_h + b_h)
    return ((h_inf - h)/tau_h)

def dmdt(m, v):
    a_m = (0.32*(v+54))/(1 - np.exp(-(v+54)/4))
    b_m = (0.28*(v+27))/(np.exp((v+27)/5)-1)
    m_inf = a_m/(a_m + b_m)
    tau_m = 1/(a_m + b_m)
    return ((m_inf - m)/tau_m)


def dndt(n, v):
    a_n = (0.032*(v+52))/(1 - np.exp(-(v+52)/5))
    b_n = 0.5*np.exp(-(v+57)/40)
    n_inf = a_n/(a_n + b_n)
    tau_n = 1/(a_n + b_n)
    return ((n_inf - n)/tau_n)


def dvdt(v, h, m, n, i_syn, i_ext): 
    Cm  = 100     
    gK = 8000 
    gL = 10 
    gNa = 10000   
    V_k = -100.0 
    V_l = -66  # leak potential (mV)
    V_na = 50.0    # Na reversal potential (mV)
    ina = gNa * h * (m**3) * (v - V_na)
    ik = gK * (n**4) * (v - V_k)
    il = gL * (v - V_l)
    return (- ina - ik - il - i_syn + i_ext )/ Cm

AMPA = {
        'g' : 7, # mS/cm2 chosen so that SC elicits more than 1 spike.
        'E' : 0,    #mV
        'a' : 1.1,   # 
        'b' : 0.19   #
        }

GABA_A = {
        #'g' : 1,
        'g':14, #maybe decrease this.
        'E' : -80,
        'a' : 5,#msec-1 mM-1. ?
        'b' : 0.18
        }

def drdt(a,b,r,NeuroTrans):
    f = (NeuroTrans * (1-r)*a - r*b)
    return  f

