clear all;
m=1; %mode number
fn='';
Nx=100; %discretize
Ny=100;

% all lengths in microns
 rwg.LAMBDA=0.6330;
 rwg.EC=1.00; % cover permittivity
 rwg.ES=3.80; % substrate permittivity
 rwg.ER=5.80; % rib permittivity
 % computational window, [xlo,xhi,ylo,yhi]
 rwg.CW=[0.0,2.5,0.0,4.0];
  % rib, [xlo,xhi,ylo,yhi]
rwg.RB=[1.25,1.75,1.5,2.5];
d=digitize(rwg,Nx,Ny);
h=helmholtz(d);
s=solve(h);
visualize(s,m,fn);