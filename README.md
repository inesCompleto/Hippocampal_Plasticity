# Hippocampal_Plasticity

The code files in this repository reproduce the modeling results in 'Recurrent cholinergic inputs induce local hippocampal plasticity through feedforward disinhibition'.

The files in the 'Cholinergic co-pairing' folder are used to study the cholinergic mechanisms by which activation of a7 nAChRs on OLM cells, paired with glutamatergic stimulation of fast-spiking interneurons and CA1 pyr, potentiates the SC-CA1 synapse. 
- ‘fast_neurons.py’, ‘OLM.py’ and ‘E.py' define the functions and parameters describing the dynamics of a fast-spiking interneuron, OLM interneuron, and dendritic compartment of a CA1 pyr cell, respectively
- 'ode.py' contains the parameters of the synaptic plasticity model, and a 'body_ode' function where the loop calculating each variable for each time step is defined, i.e. the ODE's described in the Methods section.
- 'run.py' runs the simulation for a total time of tmax=40 minutes. Everytime the user runs the code, it calculates the all the variables, and it saves the resultant EPSC and respective time step to a file named ‘EPSC_v#.txt’, where v# indicates the number of simulation the user is running.
- 'plots.py' imports all the 'EPSC_v#.txt' that the user specifies, and calculates the average and standard deviation of the EPSC of the simulations considered (for example, see Figure 2). As a practical example, we have uploaded 3 files with simulation results for a 5 minutes (EPSC_5mints_v1.txt, EPSC_5mints_v2.txt, EPSC_5mints_v3.txt), and 8 minutes co-pairing period (EPSC_8mints_v1.txt, EPSC_8mints_v2.txt, EPSC_8mints_v3.txt). Run 'plots.py' to see the corresponding graphs.
- 'run_timing' simulates changes in the co-pairing period of single pulses and respective plasticity induction (for no-noise simulations comment line ??? instead of ??? in 'ode.py'). See 'ACh_pulses.py' for examples on how to simulate paired periodic pulses (such as doublets).

The files in the 'Disinhibition co-pairing' folder are used to study the disinhibitory mechanisms by which disinhibition of the CA1 pyr cell dendritic compartment potentiates the SC-CA1 synapse when paired with SC stimulation.
- 'short_long_DisPeriod.py' runs a simulation for a short (4 mints) and long (8mints) disinhibition period. It calculates all the variables for the short and long disinhibition period and it saves the corresponding EPSC in a file named 'EPSC_8mint_v#.txt' and 'EPSC_4mint_v#.txt'. Comment line ??? instead of ?? for non-noisy simulations.
- 'plots.py' imports all the 'EPSC_8mint_v#.txt' and 'EPSC_4mint_v#.txt' that the user specifies, and calculates the average and standard deviation of the EPSC of the simulations considered.
- 'areas_ratio.py' calculates the ratio between the integral of calcium when its concentration is above the potentiation and depression onset, and the corresponding changes in AMPAR maximal conductance being induced (for an initial conductance of 8.83 nS).
