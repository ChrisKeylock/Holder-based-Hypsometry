# Holder-based-Hypsometry
Matlab code for the analysis described in Keylock et al. Water Resour. Res. 2019

A detailed description of the code is given in the commented lines at the top of the file.

Given a set of elevations (z) and corresponding Holder exponents (alpha), the code extracts information on the elevations conditioned on the Holder exponents, and the Holder exponents conditioned on the elevations. The philosophy follows that made popular in the classic paper by Arthur Strahler:  
Strahler, A. N. (1952), Hypsometric (area-altitude) analysis of erosional topography, Geol. Soc. Am. Bull., 63, 1117-1142.

In order to obtain the Holder exponents corresponding to the elevations, it is recommended that you use Fraclab:
https://project.inria.fr/fraclab/
