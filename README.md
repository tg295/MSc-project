# MSc-project
Computational models of peripheral nerve fibres (specifically axons) for conduction block - constructed in the NEURON simulation environment. 

# Final Models 

Final models folder contains the three different models constructed in the project:  

H-H model is a basic unmyelinated Hodgkin-Huxley model axon, MRG model is a myelinated axon based off of the McIntyre-Richardson-Grill (MRG) model, and the C-fibre model is an unmyelinated fibre geometrically similar to the H-H model but with added ion channels to closely represent that of C-fibres.

In both the H-H model and C-fibre model folders the axon model itself can be found in the files titled axon.hoc. The simulation protocols for stimulation are then deliniated in the init___.hoc files, including separate files for sinusoidal and square wave stimulation, as well as the GUI. In the MRG model however the simulation protocols, GUI and axon model are both contained within the init.hoc files. 

.mod contain the differential equations for ion channels and stimulation waveforms as well as additional modifications for extracellular stimulation of the axon. 

# Data Analysis 

Data analysis folder contains the scripts used in Matlab to analyse the exported results from NEURON simulations. For each model there are scripts for both individual simulations conducted using the GUI, and scripts to analyse batch runs of various stimulation protocols. 

# Thesis 

Contains the pdf document for my final report.

