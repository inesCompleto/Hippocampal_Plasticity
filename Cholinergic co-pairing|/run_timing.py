
import numpy as np
import matplotlib.pyplot as plt

import OLM
import fast_neurons as fs
import E

from ode import body_new, noise_generation

# Define time step and total time of simulation
tmin = -10
dt = 0.02

# NEUROTRANSMITTERS
Glu_dur = 5
ACh_dur = 5

init_max = 250
init_min=-30  # init and j: time difference between ACh and Glu

steps_init_time = ((init_max - init_min)/3) // 1 # number of steps
#steps_init_time = ((init_max - init_min)/10) // 1 # number of steps
#steps_init_time = ((init_max - init_min)/20) // 1
init_time = np.linspace(init_min, init_max, int(steps_init_time)) # time vector


file = open('DgDt_v9.txt','w')


delta_g = []
delta_t = []

for j in init_time:
    
    tmax=960 if j<0 else (j+960)
    
    steps = ((tmax - tmin)/dt) // 1 # number of steps
    T = np.linspace(tmin, tmax, int(steps)) # time vector

    # Initialize variables
    # Fast-spiking interneuron
    V_I = -64
    m_I=0.0
    h_I=0.0
    n_I=0.0
    r_A = 0
    r_G = 0
    GABA_I = 0

    # OLM neuron 
    V_O = -60
    Ca_a7 = 0
    Ca_is = 0.44e-3
    GABA_O = 0
    h_O = 0.0
    n_O = 0.0
    m_O = 0.0
    p = 0
    hf = 0
    hs = 0
    ra7 = 0
    ra_O=0

    # E-cell
    ra_E = 0                            #AMPAR opening gate
    rn_E = 0                            #NMDAR opening gate
    rg_E = 0                            #(GABA_A)R opening gate              
    Ca = 0             #intracellular Ca concentration
    V_Ed = -67                            #membrane potential of CA1 pyr cell            
    g = E.AMPA['g']

    # Neurotransmitter
    ACh_init = 900 # start neurotransmitter admistration at 900 msec. By doing so we make sure all variables are at equilibrium.
    Glu_init = j + 900

    x0 = (V_O, m_O, h_O, n_O, p, hf, hs, ra7, Ca_a7, Ca_is, GABA_O, 
          V_I, r_G, r_A, m_I, h_I, n_I, 
          V_Ed, ra_E, rn_E, rg_E, Ca, g, 0., 0.)

    for i,value in enumerate(T):
        
        # NEUROTRANSMITTERS
        ACh = 1 if ( ACh_init) < T[i] < (ACh_init + ACh_dur) else 0       
        Glu = 1 if ( Glu_init) < T[i] < (Glu_init + Glu_dur) else 0

        noise_O, noise_I, noise_Ed = noise_generation(T, dt) 
        
        # group arguments to be sent to integrate function
        arg = dt, ACh, Glu, noise_O[i], noise_I[i], noise_Ed[i]
    
        # integrate ODE's
        dynamics = body_new(x0, arg)
        x0 = dynamics
                
    delta_g.append(dynamics[-3] - E.AMPA['g'])
    delta_t.append(j)
    file.write("%s %s\n" % (delta_t[-1], delta_g[-1]))
    print(j)

file.close()

plt.figure(figsize=(10,7))
plt.plot(delta_t, delta_g, color='k')
plt.ylim(-0.4, 0.9)
plt.axhline(y=0, color='dimgray', linestyle='--', linewidth=1.5)
plt.axvline(x=-19.9, color='dimgray', linestyle='--', linewidth=1.5)
plt.axvline(x=10.438, color='dimgray', linestyle='--', linewidth=1.5)
plt.axvline(x=131.1, color='dimgray', linestyle='--', linewidth=1.5)
plt.axvline(x=177.4, color='dimgray', linestyle='--', linewidth=1.5)
plt.ylabel(' $\Delta \overline{g}_{AMPA}$ (mS/$cm^2$)', fontsize=25)
plt.xlabel('Time difference, $\Delta t$ (msec)', fontsize=25)
plt.xticks(fontsize=20) 
plt.yticks(np.arange(-0.3, 1, 0.2),fontsize=20) 
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.tight_layout()

plt.show()
        
