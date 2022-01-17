import numpy as np
import sk_dsp_comm.sigsys as ss

import OLM
import fast_neurons as fs
import E

# DEFINE FUNCTIONS
# Calcium dynamics:
def dCadt(Ca, i):
    alpha = 0.1
    alpha2 = 0.006
    tca = 12
    return (-alpha * alpha2 * i - Ca / tca)


def a7_dCadt(k, Ca_is, Ca, i):
    alpha = 0.05
    alpha2 = 0.0021e-3  # for simulations with noisy background current
    tca = 12
    return (-alpha * alpha2 * i + k ** 3 * (Ca_is - Ca) - Ca / tca)


# Synaptic plasticity:
def sig(x, b):
    return (np.exp(b * x) / (1 + np.exp(b * x)))

def dWdt(Ca, W):
    tau = P1 / (P2 + Ca ** P3) + P4
    learn_rate = 1 / tau

    Omega = (maxPot * sig(Ca - theta_p, tau_p) - maxDep * sig(Ca - theta_d, tau_d))

    return learn_rate * (Omega - r * (W - g0))


# periodic neurotransmitter pulse
def periodic_rect(t,init,tau,T):
    """
    The period is set by T and tau is the pulse width
    """
    x = np.zeros(len(t))
    x += ss.rect(np.mod(t - init,T),tau*2)
    return x


# DEFINE PARAMETERS
theta_d = 0.31 
tau_d = 900
theta_p = 0.34
tau_p = 900
maxPot = 0.0675 
maxDep = 0.0375


r = 0.0040
g0 = E.AMPA['g']

P1 = 1.5e-6
P2 = P1 * 1e-4
P3 = 13
P4 = 1


Tmax = 1
Vp = 2  # mV
kp = 5  # mV


def noise_generation(T, dt):
    # White noise Ed
    mean = 0
    std = 0.3
    num_samples = len(T)
    rand_Ed = np.random.normal(mean, std, size=num_samples)
    noise_Ed = dt ** 0.5 * rand_Ed

    # White noise I-cell
    mean = 0
    std = 0.1
    rand_I = np.random.normal(mean, std, size=num_samples)
    noise_I = dt ** 0.5 * rand_I

    # White noise O-cell
    mean = 0
    std = 1.1
    rand_O = np.random.normal(mean, std, size=num_samples)
    noise_O = dt ** 0.5 * rand_O
    
    return noise_O, noise_I, noise_Ed



def body_new(state, arg):
    '''
    This is a body of for-loop in jax (for fast calculations we have to precompile this part)
    Make a step of dynamics and return a tuple: (last state, history)
    State here is a tuple that consist of values of all DE-variables at the iteration step 't'. 
    '''

    dt, ACh, Glu, noise_O, noise_I, noise_Ed = arg # expand arguments; t is the index of iteration
    
    V_O, m_O, h_O, n_O, p, hf, hs, ra7, Ca_a7, Ca_is, GABA_O, V_I, r_G, r_A, m_I, h_I, n_I, V_Ed, ra_E, rn_E, rg_E, Ca, g, _, _ = state # expand current state

    
    # OLM
    m_O_new = m_O + dt * OLM.dmdt(m_O, V_O)
    h_O_new = h_O + dt * OLM.dhdt(h_O, V_O)
    n_O_new = n_O + dt * OLM.dndt(n_O, V_O)

    p_new = p + dt * OLM.dpdt(p, V_O)
    hf_new = hf + dt * OLM.dhfdt(hf, V_O)
    hs_new = hs + dt * OLM.dhsdt(hs, V_O)

    # Synaptic currents
    ra7_new = ra7 + dt * OLM.dr_a7dt(ra7, ACh)
    ia7 = 3 * ra7 * (V_O - OLM.nAChR['E'])

    # Neuron dynamics
    V_O_new = V_O + dt * OLM.dvdt(V_O, h_O, m_O, n_O, p, hf, hs, ia7, -260)  

    # CICR mechanism
    k = Ca_a7 / (Ca_a7 + 0.2e-3)
    Ca_a7_new = Ca_a7 + dt * a7_dCadt(k, Ca_is, Ca_a7, ia7)
    Ca_is_new = Ca_is + dt * (-k ** 3 * (Ca_is - Ca_a7) - (Ca_is - 0.44e-3) / 10)
    
  #  if V_O > 0:
  #      GABA_O = Tmax / (1 + np.exp(-(V_O - Vp) / kp))
  #  else:
    GABA_O = Tmax / (1 + np.exp(-(Ca_a7 - 0.00004) / 0.000001))

#     Fast-spiking

    # synaptic currents
    r_G_new = r_G + dt * fs.drdt(fs.GABA_A['a'], fs.GABA_A['b'], r_G, GABA_O)
    i_G = fs.GABA_A['g'] * r_G * (V_I - fs.GABA_A['E'])

    r_A_new = r_A + dt * fs.drdt(fs.AMPA['a'], fs.AMPA['b'], r_A, Glu)
    i_A = fs.AMPA['g'] * r_A * (V_I - fs.AMPA['E'])

    isyn = i_A + i_G

    # spiking currents
    m_I_new = m_I + dt * fs.dmdt(m_I, V_I)
    h_I_new = h_I + dt * fs.dhdt(h_I, V_I)
    n_I_new = n_I + dt * fs.dndt(n_I, V_I)

    V_I_new = V_I + dt * fs.dvdt(V_I, h_I, m_I, n_I, isyn, 0) 
    GABA_I = Tmax / (1 + np.exp(-(V_I - Vp) / kp))

#     Pyramidal Cell

    # Currents
    ra_E_new = ra_E + dt * E.drdt(E.AMPA['a'], E.AMPA['b'], ra_E, Glu)
    ia_E = g * ra_E * (V_Ed - E.AMPA['E'])
    
    rn_E_new = rn_E + dt * E.drdt(E.NMDA['a'], E.NMDA['b'], rn_E, Glu)
    B = 1 / (1 + np.exp(-0.062 * V_Ed) * E.NMDA['Mg'] / E.NMDA['kd0'])
    in_E = E.NMDA['g'] * B * rn_E * (V_Ed - E.NMDA['E'])
    
    rg_E_new = rg_E + dt * E.drdt(E.GABA_A['a'], E.GABA_A['b'], rg_E, GABA_I)
    ig_E = E.GABA_A['g'] * rg_E * (V_Ed - E.GABA_A['E'])

    isyn_E = ia_E + in_E + ig_E

    il = E.gL * (V_Ed - E.V_leak)

    # Synaptic Plasticity
    Ca_new = Ca + dt * dCadt(Ca, in_E)
    g_new = g + dt * dWdt(Ca, g)

    # Membrane Potential
    V_Ed_new = V_Ed + dt * (- il - isyn_E) / 100 + noise_Ed
    
    state = (V_O_new, m_O_new, h_O_new, n_O_new, p_new, hf_new, hs_new, ra7_new, 
             Ca_a7_new, Ca_is_new, GABA_O, V_I_new, r_G_new, r_A_new, m_I_new, h_I_new, n_I_new, 
             V_Ed_new, ra_E_new, rn_E_new, rg_E_new, Ca_new, g_new, ia_E, in_E)
   
        
    return state
