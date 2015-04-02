function plasmoninterface(ec,ef,lambda0)

k0 = 2*pi/lambda0;
%beta = k0 * sqrt(ef*ec/(ef+ec));
gamma = -k0*ef/sqrt(-ef-ec);
alphac = k0*ec/sqrt(-ef-ec);
x = linspace(-0.4,1,141);
Hy = abs(exp(gamma*x).*(x<0) + exp(-alphac*x).*(x>=0));
plot(x,Hy);
vline(0);

end
