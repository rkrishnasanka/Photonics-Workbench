function rwg=helmholtz(r)
k0=2*pi/r.LAMBDA;
hx=r.x(2)-r.x(1); % x step width
hy=r.y(2)-r.y(1); % y step width
[NX,NY]=size(r.eps'); % size of computational window
NV=NX*NY; % number of variables
prm=reshape(r.eps',NV,1); % permittivity
one=ones(NV,1);
md=-2*one*(1/hx^2+1/hy^2)/k0^2+prm; % main diagonal
xd=one/hx^2/k0^2; % side diagonal xx differentiation
yd=one/hy^2/k0^2; % side diagonal yy differentiation
hh=spdiags([yd,xd,md,xd,yd],[-NX,-1,0,1,NX],NV,NV);
% note that there are false links, remove
for n=NX:NX:NV-NX
hh(n,n+1)=0;
hh(n+1,n)=0;
end;
rwg=r;
rwg.hh=hh;
end % function helmholtz