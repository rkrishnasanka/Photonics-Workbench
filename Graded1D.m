clear all;
%TE eigenmode calculation for a graded 1D dielectric waveguide
format long;

lambda=0.6328;
k0=2*pi/lambda;
EC=1.000; %relative permittivity of the cover region
ES=4.800; %relative permittivity of the substrate region
Edelta=0.045; %permittivity profile delta parameter
W=4.00; %speread of permittivity profile


%discretization of the real axis

xmin=-1.0;
xmax=4*W;
h=0.1; %step size
x=(xmin:h:xmax)';
dim=size(x,1);

%define the permittivity profile

prm=(x<0).*EC+(x>=0).*(ES+Edelta*exp(-(x/W).^2));

%setting up the mode operator L

next=ones(dim-1,1)/h^2/k0^2;
main=-2.0*ones(dim,1)/h^2/k0^2+prm;
L=diag(next,-1)+diag(main,0)+diag(next,1);

%calcualtion of eigenvectors/eigenvalues

[evec, eval]=eig(L);

%plot only eigenvectors with epseff>epss

eff_eps=diag(eval);
guided=evec(:,eff_eps>ES);
plot(x,guided);
xlabel('Depth from surface (microns)');
ylabel('Electric field strength');



