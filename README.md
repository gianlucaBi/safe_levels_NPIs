# Online Primal-Dual Controller for the Control of Epidemic Outbreaks
This Matlab code implements a primal-dual online controller to determine safe levels of NPIs for the state of CO, United States. The precise traffic control problem is described in [1], referenced below.

# Required Software
* Matlab 2019


# Execution
There are two implementations of the controller. The files in the folder "STATE MODEL" implements a single-region model. The files in the folder "REGIONAL MODEL" implement a multi-region model fitted to the Local Public Health Agencies (LPHA) regions in Colorado.
In both folders, the file "main.m" is the main Matlab script.


# Documentation and References
[1] G. Bianchin, E. Dall'Anese, J. I. Poveda and A. Buchwald, “When can we safely return to normal? A novel method for identifying safe levels of NPIs in the context of COVID-19 vaccinations,” medRxiv preprint, apr. 2021, (medRxiv:2021.255350)


