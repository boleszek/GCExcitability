Code for generating Figures 3, 4, & 5 of the following paper:

Osinski & Kay (2015), Granule cell excitability regulates gamma and beta oscillations in a model of the olfactory bulb dendrodendritic microcircuit, (to be submitted)

Code was originally designed in MATLAB R2014b

Implementation is by Bolesław Osiński , to whom questions should be addressed
boleszek@uchicago.edu

——————————————————————————————————————————————
The zipped folder includes the following .m files:

InitNetwork_GCE.m 	   -  This function initializes network parameters
OB_network_GCE.m	   -  This is the numerical simulation of the Mitral-Granule network
ILFP_GCE .m		   -  This is a wrapper function for OB_netwrok_GCE.m. All simulation products as well as the simulated LFP are generated by this function
ParamSweep_GCE.m	   -  This is a wrapper function around ILFP_GCE.m. It evaluates the model for every combination of the two parameter inputs
Rasterplot.m    	   -  Generates raster plot
——————————————————————————————————————————————


Usage:
1. Unzip GCE_ModelDB.zip
2. Open MATLAB. Set MATLAB path to unzipped folder’s directory
3. Open Fig3, Fig4, or Fig5 directory and run the associated script (i.e. Fig3script.m) for producing the figures. Each folder includes it’s own parameter file (i.e. OB_params_GCE_Fig3.txt)

Note: Fig3script.m runs quickly, Fig5script.m runs for ~ 4min, and Fig4script.m runs for ~ 3.5hr

