# DiffusionFEM

Internal finite element code using Fenics. Disclaimer: I am still in the learning stage and code here may be of poor quality or just straight up incorrect. Note that many of the parameter values in these scripts (e.g. diffusivities) are approximate values based on literature data and personal experiments; however, "physically", they don't mean very much.

## Directory and File information

/ExperimentalRecreation  
Everything in this directory is meant to recreate (to varying levels of accuracy) crystal dissolution experiments that were done at UC Santa Barbara

/ExperimentalRecreation/DiffusionDirANDNeu.py  
Simple 2D diffusion solution that has a variable concentration boundary condition on one end, and no flux out of all other boundaries. 

/ExperimentalRecreation/1D  
This directory contains 1D simulations of experiments

/ExperimentalRecreation/1D/Fitting.py  
This is a simple test of numerical fitting using a FEM Fenics backend for optimizing the parameters. Example data is given based on actual experiments. This assumes that the diffusion is 1D in nature (verified experimentally)


/ExperimentalRecreation/1D/1DExperimentalRecreation.py  
Simple code meant to recreate experiments in 1D. Takes in user input for diffusivity and interface concentration


---

/NonlinearD  
Everything in this directory has code for a diffusivity that is not a function of concentration. Typically this is assumed to have some exponential variance (which is what was assumed here)

/NonlinearD/DiffFConc.py  
This is the main function for 2D diffusion with a concentration dependent diffusivity

/NonlinearD/1D/code.py  
This is code for 1D diffusion with a concentration dependent diffusivity

/NonlinearD/1D/curvefitting.py  
This code fits some given experimental data assuming that the diffusivity is a function of concetration. The fitting method again uses Fenics FEM as the way to optimize the fitting parameters, which are Do (diffusivity exponential prefactor), a (the exponential variable), and Cs (the interface concentration). 


---


/SimpleTests  
Everything in this directory is simple scripts that acts mostly as a testing ground for me. 

/SimpleTests/1Ddiffusion.py  
The first implimentation of 1D diffusion. This code served as the back bones for more complex scripts later

/SimpleTests/KillAtConcentration.py  
This code is for 2D non-linear diffusion. The user enters in the diffusivity (or fixes it in the code) and the time step. The code will iterate with a while loop until a "kill" condition is meant. This code is used to find how long it takes to saturate a model intercolumnar gap (2D - assuming no melt movement). Thus, the kill condition is typically something like 90% of the final concentration at the edges. The code then outputs the time needed. Do not use this to generate gifs of the process but this is good for useful data.

/SimpleTests/MyDiffusion.py  
Very early code - simple 2D diffusion with dissolution from all edges. Later used to try and learn outputting directly to CSV files from python.

