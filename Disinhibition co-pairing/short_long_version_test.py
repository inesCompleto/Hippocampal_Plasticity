import numpy as np
from math import exp

# ----------  PARAMETERS  ---------- #
#Synaptic Currents
AMPA = {
        'g' : 4, #maximal conductance in nS
        'E' : 0,   #reversal potential in mV
        'a' : 1.1,  #opening rate in mM-1 msec-1
        'b' : 0.19  #closing rate in msec-1
        }

NMDA = {
        'g' : 25,
        'Mg' : 1,    #Mg concentration in mM
        'kd0': 3.57, #
        'E' : 0,
        'a' : 0.072,  
        'b' : 6.6e-3 ,
        }

GABA_A = {
        'g' : 7,
        'E' : -80,
        'a' : 5,
        'b' : 0.18
        }


gL = 1  #leak conductance nS
V_leak = -68 #leak reversal potential

#Calcium Dynamics
alpha = 0.1 #fraction of NMDA current composed of Ca ions
alpha2 = 0.045
tca = 12     #decay constant in msec

#Synaptic Plasticity
theta_d = 0.31    #depression threshold in mM
tau_d = 900     # 
theta_p = 0.34  #potentiation threshold in mM
tau_p = 900
maxDep = 0.0375 
maxPot=0.0699

P1 = 1.5e-6 
P2 = P1*1e-4
P3=13
P4 = 1
r=0.0040
g0 = AMPA['g'] #initial value of gA



# ----------  FUNCTIONS  ---------- #
#Neurotransmitters

Glu_par ={
        'init': 0, 
        'dur' : 1.0,   #msec
        'max' : 1,     #mM
        'k' : 0
        }

GABA_par = {
       'init': 2,
       'dur': 1.0,   #msec
       'max': 1,     #mM
       'k': 0
       }

period = 60e3 


#Calcium Dynamics
def dCadt(Ca, i ): 
    return (-alpha2*alpha * i -  Ca/tca )

#Synaptic plasticity dynamics
def dWdt(Ca, W):    
    tau = P1/(P2 + Ca**P3) + P4
    learn_rate = 1/tau
    Omega = ( maxPot*sig(Ca - theta_p, tau_p) - maxDep*sig(Ca - theta_d, tau_d) )
    
    return learn_rate*(Omega - r*(W-g0)) 

#Auxiliar function for synaptic plasticity
def sig(x, b):
    return (exp(b*x)/(1 + exp(b*x)))

#Receptor opening gate dynamics
def drdt(a,b,r,NeuroTrans):
    if NeuroTrans == 0:
        f = - r*b
    else:
        f = (NeuroTrans * (1-r)*a - r*b)
    return  f



# ----------  VARIABLES  ---------- #
tmin = 0.0
tmax = 2.7e6
dt=0.02
steps = ((tmax - tmin)/dt) // 1
T = np.linspace(tmin, tmax, int(steps)) #total time in msec
T_min = [x*1e-3/60 for x in T]     #total time in min (for plotting)

Glu = np.zeros(len(T))

#short disinhibition 
rGb_s=0
s_l=0
rGb_l=0
s_s=0
r_A_s = 0                            #AMPAR opening gate
r_N_s = 0                            #NMDAR opening gate
r_G_s = 0                            #(GABA_A)R opening gate              
g_s = np.zeros(len(T))                     #AMPAR maximal conductance
g_s[0] = AMPA['g']
Ca_s = np.zeros(len(T))             #intracellular Ca concentration
EPSC_x_s=[]                          #time of EPSC local maximum
EPSC_y_s=[]                          #EPSC amplitude (local maximum)
V_s = -67                            #membrane potential of CA1 pyr cell            

#long disinhibition 
r_A_l = 0                            #AMPAR opening gate
r_N_l = 0                            #NMDAR opening gate
r_G_l = 0                            #(GABA_A)R opening gate              
g_l = np.zeros(len(T))                      #AMPAR maximal conductance
g_l[0] = AMPA['g']
Ca_l = np.zeros(len(T))             #intracellular Ca concentration
EPSC_y_l=[]                          #EPSC amplitude (local maximum)
V_l = -67                            #membrane potential of CA1 pyr cell  


i=0
j=0

#White Noise
mean = 0
std=40
num_samples = len(T)
rand = np.random.normal(mean,std,size=num_samples)
sigma=10
noise = sigma * dt **0.5 * rand

t_dis_s=546e3
t_dis_l=726e3


#Opening files to write on
file_s = open('EPSC_5mint.txt','w')
file_l = open('EPSC_8mint.txt','w')

# ----------  INTEGRATE VARIABLES  ---------- #
while i < len(T)-1:
    
    #Neurotransmitter Concentration
    if T[i] > ( GABA_par['init'] + GABA_par['dur'] + period*GABA_par['k'] ) :
        GABA_par['k'] += 1
    if ( GABA_par['init'] + 60e3*GABA_par['k']) < T[i] < (GABA_par['init'] + GABA_par['dur'] + 60e3*GABA_par['k']) and (T[i]<246e3 or T[i]>t_dis_s) :
        GABA_s = GABA_par['max']
    else:
        GABA_s = 0
        
    if ( GABA_par['init'] + 60e3*GABA_par['k']) < T[i] < (GABA_par['init'] + GABA_par['dur'] + 60e3*GABA_par['k']) and (T[i]<246e3 or T[i]>t_dis_l) :
        GABA_l = GABA_par['max']
    else:
        GABA_l = 0
    
    if T[i] > ( Glu_par['init'] + Glu_par['dur'] + period*Glu_par['k'] ) :
        Glu_par['k'] += 1
    if ( Glu_par['init'] + period*Glu_par['k']) < T[i] < (Glu_par['init'] + Glu_par['dur'] + period*Glu_par['k']) :
        Glu[i] = Glu_par['max']
    else: 
        Glu[i] = 0
    
    
    
    #short disinhibition 
    #Currents
    r_A_s = r_A_s + dt*drdt(AMPA['a'], AMPA['b'], r_A_s, Glu[i])
    i_A_s = g_s[i] * r_A_s * (V_s - AMPA['E'])
    r_N_s = r_N_s + dt*drdt(NMDA['a'], NMDA['b'], r_N_s, Glu[i])
    B_s = 1 / (1 + np.exp(-0.062 * V_s) * NMDA['Mg']/NMDA['kd0']) 
    i_N_s = NMDA['g'] * B_s * r_N_s * (V_s - NMDA['E'])
    r_G_s = r_G_s + dt*drdt(GABA_A['a'], GABA_A['b'], r_G_s, GABA_s)
    i_G_s = GABA_A['g'] * r_G_s * (V_s - GABA_A['E'])
    isyn_s = i_A_s + i_N_s + i_G_s
    
    il_s = gL * (V_s - V_leak)
    
    #Synaptic Plasticity 
    Ca_s[i+1] = Ca_s[i] + dt * dCadt(Ca_s[i], i_N_s) 
    g_s[i+1] =  g_s[i] + dt *dWdt(Ca_s[i], g_s[i])  
    
    #Membrane Potential
    V_s = V_s + dt *  (- il_s - isyn_s)/100  
#    V_s = V_s + dt *  (- il_s - isyn_s + noise[i] )/100 
    
    #EPSC
    EPSC_s = -i_A_s - i_N_s 
    
    #long disinhibition 
    #Currents
    r_A_l = r_A_l + dt*drdt(AMPA['a'], AMPA['b'], r_A_l, float(Glu[i]))
    i_A_l = g_l[i] * r_A_l * (V_l - AMPA['E'])
    r_N_l = r_N_l + dt*drdt(NMDA['a'], NMDA['b'], r_N_l, float(Glu[i]))
    B_l = 1 / (1 + np.exp(-0.062 * V_l) * NMDA['Mg']/NMDA['kd0']) 
    i_N_l = NMDA['g'] * B_l * r_N_l * (V_l - NMDA['E'])
    r_G_l = r_G_l + dt*drdt(GABA_A['a'], GABA_A['b'], r_G_l, float(GABA_l))
    i_G_l = GABA_A['g'] * r_G_l * (V_l - GABA_A['E'])
    isyn_l = i_A_l + i_N_l + i_G_l  
    
    il_l = gL * (V_l - V_leak)
    
    #Synaptic Plasticity 
    Ca_l[i+1] = Ca_l[i] + dt * dCadt(Ca_l[i], i_N_l) 
    g_l[i+1] =  g_l[i] + dt *dWdt(Ca_l[i], g_l[i])
    
    #Membrane Potential
    V_l = V_l + dt *  (- il_l - isyn_l)/100  
#    V_l = V_l + dt *  (- il_l - isyn_l + noise[i] )/100 
    
    #EPSC
    EPSC_l = -i_A_l - i_N_l
    
    
    if Glu[i]>Glu[i-1]:
        EPSC_x_s.append(T_min[i])
        EPSC_y_s.append(EPSC_s)
        
        EPSC_y_l.append(EPSC_l)
        
        file_s.write("%s %s\n" % (EPSC_x_s[j], EPSC_y_s[j]))
        file_l.write("%s %s\n" % (EPSC_x_s[j], EPSC_y_l[j]))
        print(j)
        j+=1  

    i+=1

file_s.close()
file_l.close()

# ----------  PLOTS  ---------- #
import matplotlib.pyplot as plt
plt.figure()
plt.plot(EPSC_x_s, EPSC_y_s)
plt.plot(EPSC_x_s, EPSC_y_l)

plt.show()
