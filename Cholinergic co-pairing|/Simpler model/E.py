"""

Assuming surface area of 1e-5 cm2:
    pF, nS, pA, mV, msec.

"""



def drdt(a,b,r,NeuroTrans):
    f = (NeuroTrans * (1-r)*a - r*b)
    return  f

# synaptic currents function
AMPA = {
        'g' : 4, #maximal conductance in nS
        'E' : 0,   #reversal potential in mV
        'a' : 1.1,  #opening rate in mM-1 msec-1
        'b' : 0.19  #closing rate in msec-1
        }

NMDA = {
        'g' : 25,
   #'g' : 5,
        'Mg' : 1,    #Mg concentration in mM
        'kd0': 3.57, #
        'E' : 0,
        'a' : 0.072,  
        'b' : 6.6e-3 ,
        }

GABA_A = {
        #'g' : 8,
        'g' : 7,
        'E' : -80,
        'a' : 5,
        'b' : 0.18
        }


gL = 1  #leak conductance nS
V_leak = -68 #leak reversal potential
