# Hippocampal_Plasticity

The code files in this repository reproduce the modeling results in 'Recurrent cholinergic inputs induce local hippocampal plasticity through feedforward disinhibition'.

The files in the 'Cholinergic co-pairing' folder are used to study the cholinergic mechanisms by which activation of a7 nAChRs on OLM cells, paired with glutamatergic stimulation of fast-spiking interneurons and CA1 pyr, potentiates the SC-CA1 synapse. 
- ‘fast_neurons.py’, ‘OLM.py’ and ‘E.py' define the functions and parameters describing the dynamics of a fast-spiking interneuron, OLM interneuron, and dendritic compartment of a CA1 pyr cell, respectively
- 'ode.py' contains the parameters of the synaptic plasticity model, and a 'body_ode' function where the loop calculating each variable for each time step is defined, i.e. the ODE's described in the Methods section. For noisy simulations, consider line ????; for non-noisy simulations, consider line ??? instead.
- 'run.py' runs the simulation for a total time of tmax=40 minutes. Everytime the user runs the code, it calculates the all the variables, and it saves the resultant EPSC and respective time step to a file named ‘EPSC_v#.txt’, where v# indicates the number of simulation the user is running.
- 'plots.py' imports all the 'EPSC_v#.txt' that the user specifies, and calculates the average and standard deviation of the EPSC of the simulations considered (for example, see Figure 2).
- 'run_timing' reproduces Figure 3.

The files in the 'Disinhibition co-pairing' folder are used to study the disinhibitory mechanisms by which disinhibition of the CA1 pyr cell dendritic compartment potentiates the SC-CA1 synapse when paired with SC stimulation.
