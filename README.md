# DiffusionFEM

Internal finite element code using Fenics. Disclaimer: I really don't know what I'm doing :)

## DiffusionDirANDNeu.py

This code is meant to capture a constant diffusivity (basic diffusion PDE) of TBC materials in model silicate melts. The boundary conditions are a combination of Dirichlet (time-evolving) in a small region of the boundary, and are zero flux out of the boundaries (Neumann) everywhere else; these boundary coditions are meant to capture experimental conditions.

## DiffFConc.py

This code is meant to capture a non-constant diffusivity (diffusivity that is a function of concentration). The boundary conditions are setup to be a 2D recreation of my experiments. Here, we are modeling the diffusion of an element (e.g. ZrO2) in model silicate melts.

## curvefitting.py

This code is meant to use the same elements of DiffFConc.py in a 1-D application for use in curve fitting experimental data. Example data is provided in the code (xdata, ydata). The curve fitting is done using scipy and function generation is done with Fenics.
