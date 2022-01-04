"""

Assuming surface area of 1e-5 cm2:
    pF, nS, pA, mV, msec.

"""


import numpy as np

def dhdt(h, v):
    a_h = 0.07*np.exp(-(v+37)/20)
    b_h = 1/(np.exp(-0.1*(v+7))+1)
    h_inf = a_h/(a_h + b_h)
    tau_h = 1/(a_h + b_h)
    return ((h_inf - h)/tau_h)

def dmdt(m, v):
    a_m = (-0.1*(v+23))/(np.exp(-0.1*(v+23)) - 1)
    b_m = 4*np.exp(-(v+48)/18)
    m_inf = a_m/(a_m + b_m)
    tau_m = 1/(a_m + b_m)
    return ((m_inf - m)/tau_m)


def dndt(n, v):
    a_n = (-0.01*(v+27))/(np.exp(-0.1*(v+27))-1)
    b_n = 0.125*np.exp(-(v+37)/80)
    n_inf = a_n/(a_n + b_n)
    tau_n = 1/(a_n + b_n)
    return ((n_inf - n)/tau_n)


def dpdt(p, v):
    a_p = 1/(0.15*(1+np.exp(-(v+38)/6.5)))
    b_p = np.exp(-(v+38)/6.5)/(0.15*(1+np.exp(-(v+38)/6.5)))
    p_inf = a_p/(a_p + b_p)
    tau_p = 1/(a_p + b_p)
    return ((p_inf - p)/tau_p)


def dhfdt(hf, v):
    hf_inf = 1/(1+np.exp((v+79.2)/9.78))
    tau_hf = 1 + 0.51/(np.exp((v-1.7)/10) + np.exp(-(v+340)/52))
    return ((hf_inf - hf)/tau_hf)

def dhsdt(hs, v):
    hs_inf = 1/(1+np.exp((v+2.83)/15.9))**58
    tau_hs = 5.6/(np.exp((v-1.7)/14) + np.exp(-(v+260)/43)) + 1
    return ((hs_inf - hs)/tau_hs)


def dvdt(v, h, m, n, p, hf, hs, i_syn, i_ext): 
    Cm  = 100  
    gp = 50
    gh = 145
    V_h = -20 
    gK = 1100 
    gL = 50
    gNa = 5200  
    V_k = -90 
  #  V_l = -65  # leak potential (mV)
    V_l = -70    # different from Rotstein to have resting potential at -60, as reported by Leão.
    V_na = 55    # Na reversal potential (mV)
    
    ina = gNa * h * (m**3) * (v - V_na)
    ik = gK * (n**4) * (v - V_k)
    il = gL * (v - V_l)
    ip = gp * p * (v - V_na)
    ih = gh * (0.65*hf + 0.35*hs)*(v - V_h) 

    return (- ina - ik - il - ih - ip - i_syn + i_ext )/ Cm



# synaptic currents function
nAChR = {
        'g' : 3, #chosen to have activation of a7 nAChR with a pulse of ACh evoke a current with 35pA (Leão et al.)
        'E' : 0 
       }

def dr_a7dt(r, NeuroTrans):
    n = 1.73 # Hill's coefficient of activation
    EC50 = 80e-3 #mM, half-maximum concentration of activation
    tau = 5 # msec, activation time cte
    a_inf = NeuroTrans**n / (EC50**n + NeuroTrans**n)
    return ((a_inf - r)/tau)


def dr_s7dt(r, NeuroTrans):
    n = 2
    IC50 = 1.3e-3 #mM, half-maximum concentration of activation
    tau_0 = 50 # msec, activation time cte
    tau_max = 120e3
    K_tau = 1.73e-3
    eta=1
    a_inf = IC50**n / (IC50**n + (eta*NeuroTrans)**n)
    tau = tau_0 + tau_max*K_tau**n/(K_tau**n + (eta*NeuroTrans)**n)
    a = a_inf / tau
    b = 1/tau
    return (a * (1 - r) - b * r)


AMPA = {
    'g' : 1,
        'E' : 0,   #reversal potential in mV
        'a' : 1.1,  #opening rate in mM-1 msec-1
        'b' : 0.19  #closing rate in msec-1
        }
