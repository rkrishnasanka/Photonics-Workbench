function result=f_TE(z)
%creates function needed to calculate the propagation constant with the
%transfer method
% result - expression used in search for the propagation constant
% z - actual value of the propagation constant

global n_s
global n_c
global n_layer
global d_layer
global k_0

zz=z*k_0;
NumLayers=length(d_layer);

%creation of substrate and cladding
gamma_sub=sqrt(zz^2-(k_0*n_s)^2);
gamma_clad=sqrt(zz^2-(k_0*n_c)^2);

%creation of kappa for internal layers
kappa=sqrt(k_0^2*n_layer.^2-zz.^2);
temp=kappa.*d_layer;

%construction of transfewr matrix for first layer

cc=cos(temp);
ss=sin(temp);
m(1,1)=cc(1);
m(1,2)=-1j*ss(1)/kappa(1);
m(2,1)=-1j*kappa(1)*ss(1);
m(2,2)=cc(1);

%construction of transfer matrix for the remaining layers

for i=2:NumLayers
    mt(1,1)=cc(i);
    mt(1,2)=-1j*ss(i)/kappa(i);
    mt(2,1)=-1j*kappa(i)*ss(i);
    mt(2,2)=cc(i);
    m=mt*m;
end

result=1j*(gamma_clad*m(1,1)+gamma_sub*m(2,2))+m(2,1)-gamma_sub*gamma_clad*m(1,2);

