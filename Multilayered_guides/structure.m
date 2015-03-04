%contains data for a general waveguide

%global variables to be transferred to function f_TE.m

global n_c      %ref index cladding
global n_layer  %ref index of internal layers
global n_s      %refr index of substrate
global d_c      %thickness of cladding (microns)
global d_layer  %thickness of internal layers (microns)
global d_s      %thickness of substrate (microns)
global k_0       %wavenumber
global NumberMesh %number of mesh points in each layers including substrate and cladding


n_c=1.0;
n_layer=[1.66, 1.60, 1.53, 1.66]; %internal layers
n_s=1.5;
d_c=0.5;
d_layer=[0.5 0.5 0.5 0.5]; %thicknesses of internal layers
d_s=1.0;
NumberMesh=[10 10 10 10 10 10]; %size equal to internal layer number plus substrate and cladding!
lambda=0.6328; %wavelength in mocrons
k_0=2*pi/lambda;