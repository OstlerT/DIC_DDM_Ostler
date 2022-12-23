## Invariance of Differential Dynamic Microscopy analysis in Differential Interference Contrast images
 Code to generate data and figures from the Ostler et al DIC/DDM paper


# What is this code?
This code provides all of the supporting simulations and analysis for the named paper. The code for DDM analysis presented here been adapted from code made available by [Germain et al](https://doi.org/10.1119/1.4939516), to whom the authors are most grateful.

We also present code which simulates particle motion and generates corresponding images, for DDM to be applied to. Although this code is written by Ostler et al, the simulation parameters, and choice of Point Scattering Function for the appearance of the particles, is taken from work by [Bayles et al](DOI	https://doi.org/10.1039/C5SM02576A), in order to ensure our simulations are comparable to those in literature.

The main file in this repository is DIC_DDM, which runs the user through a loop of generating data, and allows users to generate some figures from the paper. The file 'Simulation_Functions' should be added to path, as it contains prerequisite functions to run the original simulations and analysis.
#Prerequisites
This code has been run in Matlab 2022b, and requires the Signal Processing Toolbox. Users without access to this can edit line 10 in DDMAlgorithm.mat to replace the Blackman-Harris filter with an equivalent alternative.


Colloidal dataset found at https://doi.org/10.6084/m9.figshare.21777140.