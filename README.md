# EyeMouseInteraction
Calculating interaction between the eye and mouse trajectories from an eye-tracking experiment in HCI

This code is free and open and provided as supplementary information to the following article:

Demšar U and Çöltekin A, 2017, Quantifying Gaze and Mouse Interactions on Spatial Visual Interfaces with a New Movement Analytics Methodology. PloS One 12(8): e0181818.
https://doi.org/10.1371/journal.pone.0181818

The code is provided as is and can be changed, but we would appreciate if you would please, cite the paper 
if you use the code as well as send us a message - we will be delighted if anyone is interested 
in our work.

Code uses and implements methods for stacked space-time densities and Brownian bridges. 

For stacked space-time density algorithm see:

Demsar U, Buchin K, van Loon EE and Shamoun-Baranes J, 2015, Stacked space-time densities: a geovisualisation approach to explore 
dynamics of space use over time. GeoInformatica, 19(1):85-115. DOI 10.1007/s10707-014-0207-5

For algorithm for Brownian bridge decay see:
Horne JS, Garton EO, Krone SM, Lewis JS, 2007, Analyzing animal movements using Brownian bridges, Ecology 88(9):2354–2363

# How to run the code:

The main file is DynamicInteractionEyeAndMouse.R, which specifies values for input files,
parameters and outputs. The file runs on example eye/mouse data from one task and two participants, 
which are provided as part of this zip archive. Parameters are set for foveal interaction.

For use on other data and with other input parameters, the top section of the file DynamicInteractionEyeAndMouse.R file ('Parameters to be set/changed by the user') 
should be appropriately edited.
