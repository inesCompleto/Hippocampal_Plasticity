# Hippocampal_Plasticity

The code files in this repository reproduce the modeling results in 'Recurrent cholinergic inputs induce local hippocampal plasticity through feedforward disinhibition'.

The files in the 'Cholinergic co-pairing' folder are used to study the cholinergic mechanisms by which activation of a7 nAChRs on OLM cells, paired with glutamatergic stimulation of fast-spiking interneurons and CA1 pyr, potentiates the SC-CA1 synapse. 
- ‘fast_neurons.py’, ‘OLM.py’ and ‘E.py' define the functions and parameters describing the dynamics of a fast-spiking interneuron, OLM interneuron, and dendritic compartment of a CA1 pyr cell, respectively
- 'ode.py' contains the parameters of the synaptic plasticity model, and a 'body_ode' function where the loop calculating each variable for each time step is defined, i.e. the ODE's described in the Methods section.
- 'run.py' runs the simulation for a total time of tmax=40 minutes. Everytime the user runs the code, it calculates the all the variables, and it saves the resultant EPSC and respective time step to a file named ‘EPSC_v#.txt’, where v# indicates the number of simulation the user is running.
- 'plots.py' imports all the 'EPSC_v#.txt' that the user specifies, and calculates the average and standard deviation of the EPSC of the simulations considered (for example, see Figure 2). As a practical example, we have uploaded 3 files with simulation results for a 5 minutes (EPSC_5mints_v1.txt, EPSC_5mints_v2.txt, EPSC_5mints_v3.txt), and 8 minutes co-pairing period (EPSC_8mints_v1.txt, EPSC_8mints_v2.txt, EPSC_8mints_v3.txt). Run 'plots.py' to see the corresponding graphs.
- 'run_timing' reproduces Figure 3C (for no-noise simulations comment line ??? instead of ??? in 'ode.py').
- 'run_timing_pulses' reproduces Figure 3D (for no-noise simulations comment line ??? instead of ??? in 'ode.py').

The files in the 'Disinhibition co-pairing' folder are used to study the disinhibitory mechanisms by which disinhibition of the CA1 pyr cell dendritic compartment potentiates the SC-CA1 synapse when paired with SC stimulation.
