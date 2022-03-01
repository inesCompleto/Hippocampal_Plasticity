import numpy as np
from math import exp
from scipy import integrate


# ----------  PARAMETERS  ---------- #
#Synaptic Currents
AMPA = {
        'gS' : 8.83, #maximal conductance in mS/mm2
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
       # 'g' : 0.2,
        #'g' : 10,
        'g' : 7,
        'E' : -80,
        'a' : 5,
        'b' : 0.18
        }


gL = 1  #leak conductance nS
V_leak = -68 #leak reversal potential

# ----------  FUNCTIONS  ---------- #
def dCadt(Ca, i, iCa2 ): 
    alpha=0.1
    alpha2=0.045
    tca=12

    return (-alpha2*alpha * i - alpha2*iCa2 -  Ca/tca )

#Synaptic plasticity dynamics
def dWdt(Ca, W):
    theta_d = 0.31
    tau_d = 900 
    theta_p = 0.34
    tau_p = 900
    maxPot=0.0699
    maxDep=0.0375
    
    P1 = 1.5e-6
    P2 = P1*1e-4
    P3 = 13
    P4 = 1
    r=0.0040
#    g0 = AMPA['g']
    
    tau = P1/(P2 + Ca**P3) + P4
    learn_rate = 1/tau
    
    Omega = ( maxPot*sig(Ca - theta_p, tau_p) - maxDep*sig(Ca - theta_d, tau_d) )
    
    return learn_rate*(Omega - r*(W-4))

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
tmax=650
dt = 0.01
steps = ((tmax - tmin)/dt) // 1
T = np.linspace(tmin, tmax, int(steps)) #total time in msec



Glu = []
GABA=[]
for t in T:
    if 0<t<1:
        Glu.append(1)
    else:
        Glu.append(0)

for t in T:
    if 2<t<3:
        GABA.append(1)
    else:
        GABA.append(0)



r_A = 0                            #AMPAR opening gate
r_N = 0
r_G = 0  #(GABA_A)R opening gate 
V = -67                           #membrane potential of CA1 pyr cell            
Ca = np.zeros(len(T))
g = np.zeros(len(T))    
g[0] = AMPA['gS']

area_pot=0
area_dep=0

i=0

onset=51
tau=25
iCa=np.zeros(len(T))
Omega=np.zeros(len(T))
learn_rate=np.zeros(len(T))
weight=np.zeros(len(T))
theta_d = 0.31
tau_d = 900 
theta_p = 0.34
tau_p = 900
maxPot=0.0699
maxDep=0.0375
  
P1 = 1.5e-6
P2 = P1*1e-4
P3 = 13
P4 = 1
  

func=np.zeros(len(T))

# ----------  INTEGRATE VARIABLES  ---------- #
while i < len(T)-1:
    

    #short disinhibition 
    #Currents
    r_A = r_A + dt*drdt(AMPA['a'], AMPA['b'], r_A, Glu[i])
    i_A = (g[i] * r_A * (V - AMPA['E']))
    r_N = r_N + dt*drdt(NMDA['a'], NMDA['b'], r_N, Glu[i])
    B = 1 / (1 + np.exp(-0.062 * V) * NMDA['Mg']/NMDA['kd0']) 
    i_N = (NMDA['g'] * B * r_N * (V - NMDA['E']))
    r_G = r_G + dt*drdt(GABA_A['a'], GABA_A['b'], r_G, GABA[i])
    i_G = (GABA_A['g'] * r_G * (V - GABA_A['E']))
    isyn = i_A + i_N + i_G
    
    
    il = gL * (V - V_leak)
    
    if T[i]>onset:
        gCa = 0.0023*(T[i] - onset)/tau*exp(-(T[i] - onset)/tau)
        iCa[i] = gCa*(V - 140)
    else:
        iCa[i]=0
  
    
    #Synaptic Plasticity 
    Ca[i+1] = Ca[i] + dt * dCadt(Ca[i], i_N, 0) 
    g[i+1] =  g[i]  + dt *dWdt(Ca[i], g[i])
    
    Omega[i+1] = ( maxPot*sig(Ca[i+1] - theta_p, tau_p) - maxDep*sig(Ca[i+1] - theta_d, tau_d) )
    learn_rate[i+1] = 1/(P1/(P2 + Ca[i+1]**P3) + P4)     
    
    #Membrane Potential
    V = V + dt *  (- il - isyn)/100 
    
    # Calculate area of potentiation using trapezoidal rule
    weight[i+1] = learn_rate[i+1]
    func[i] = weight[i]*Ca[i]
    
    if Ca[i]>0.34:
        area_pot+= 0.5*(T[i+1]-T[i])*(Ca[i+1]*weight[i+1] + Ca[i]*weight[i] )
    
    if 0.31<Ca[i]<0.34:
        area_dep+= 0.5*(T[i+1]-T[i])*(Ca[i+1]*weight[i+1] + Ca[i]*weight[i] )
    
    i+=1


print(area_pot/area_dep)
print(g[-1] - g[0])  



# ----------  PLOTS  ---------- #

import matplotlib.pyplot as plt


plt.figure()
plt.plot(T,Ca, color='k')
#plt.plot(T,func)
plt.hlines(y=0.31, xmin=0, xmax=900, color='dimgray', linestyle='--')
plt.hlines(y=0.34, xmin=0, xmax=900, color='dimgray', linestyle='--')
#plt.hlines(y=0.36, xmin=0, xmax=900, color='dimgray', linestyle='--')
plt.text(301, 0.309, '$\Theta_d$', fontsize=15)
plt.text(301, 0.339, '$\Theta_p$', fontsize=15)
#plt.text(301, 0.359, '$\Theta_p^{eff}$', fontsize=15)
plt.text(200, 0.55, '$\\bar{g}_{AMPA} = 8.83\ nS$', fontsize=20)
plt.text(200, 0.50, '$\\frac{A_{pot}}{A_{dep}}=9.25$', fontsize=20)
plt.ylim(0,0.55)
plt.xlim(0,300)
plt.fill_between(T,Ca,0, where = (Ca>=0.34), color='#B96641', alpha=0.6, label='$A_{pot}$')
plt.fill_between(T,Ca,0, where = (Ca<0.34) & (Ca>0.31), color='#505D74', alpha=0.6, label='$A_{dep}$')
plt.ylabel('Calcium ($\mu$M)', fontsize=20)
plt.xlabel('time (msec)', fontsize=20, x=1,y=1)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
#plt.legend(loc=2,fontsize=20)


plt.figure()
plt.plot(Ca,g, color='k')
#plt.vlines(x=0.31, ymin=6, ymax=8, color='dimgray', linestyle='--')
#plt.vlines(x=0.34, ymin=6, ymax=8, color='dimgray', linestyle='--')
plt.ylim(6,8)
plt.annotate("", xy=(Ca[501], g[501]), xytext=(Ca[500], g[500]), arrowprops=dict(arrowstyle="->",mutation_scale=20,linewidth=1.5, color='k'))
plt.annotate("", xy=(Ca[2001], g[2001]), xytext=(Ca[2000], g[2000]), arrowprops=dict(arrowstyle="->",mutation_scale=20,linewidth=1.5, color='k'))
plt.annotate("", xy=(Ca[6001], g[6001]), xytext=(Ca[6000], g[6000]), arrowprops=dict(arrowstyle="->",mutation_scale=20,linewidth=1.5, color='k'))
plt.annotate("", xy=(Ca[15501], g[15501]), xytext=(Ca[15500], g[15500]), arrowprops=dict(arrowstyle="->",mutation_scale=20,linewidth=1.5, color='k'))


plt.show()



