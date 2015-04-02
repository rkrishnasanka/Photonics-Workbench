%%=========================================================================
% EC568 fiber modesolver
% Original by Siddharth Ramachandran; Modified by Yuhao Chen
% Yuhao Last Update: Oct. 30, 2011
% Siddharth Last Update: March 30, 2015
% NOTE: EVERYTHING (EXCEPT SELLEMEIR COEFFS) IS IN SI UNITS
%=========================================================================

clear all
close all

% NOTE: EVERYTHING (EXCEPT SELLEMEIR COEFFS) IS IN SI UNITS
c=3e8; %(m/s)

% Make physical dimensions of the fiber (radius)
    r1=10e-6;   % Core (m)
    r2=60e-6;  % Cladding (m)
    del_r=0.01e-6; % resolution of the radius of the fiber cross section
    
% Specify NA of fiber. Code automatically calculates n1 assuming n2 = n_SiO2 
    NA = 0.1;
    
% N = Number of points of the index profile over which computation is done.
% This is the # of points in the field plot that you get at the output.
    N = 300;
         
%Solves from LP(L,0) to LP(L,m) modes
    L = 0;
    m = 3; 
   
%And now we specify the wavelength range over which you want the simulation
%to run: lambdastart and lambdaend, and the wavelength spacing: lambdastep.
   lambdastart=0.7e-6; %(m)
   lambdaend=1.7e-6;   %(m)
   lambdastep=0.02e-6; %(m)
   lambda=lambdastart:lambdastep:lambdaend;
   w=2.*pi.*c./lambda; 
   
    
%% ========================================================================
% PROGRAM STARTS
%Initialize wavelength loop

% Pure silica SiO2 coefficients in sellemeier's formula (Not in SI units)
    a1=0.004679148;
    a2=0.01351206;
    a3=97.93400;
    b1=0.6961663;
    b2=0.4079426;
    b3=0.8974794;

% dont change ... building rr profile
    rr1=0:del_r:r1;
    rr2=r1+del_r:del_r:r2;
    rr=cat(2,rr1,rr2);
    
% initialise parameters space  
n_material=zeros(1,numel(lambda));
beta=zeros(numel(lambda),m);
neff=zeros(numel(lambda),m);
field=zeros(N,m,numel(lambda));
Dispersion=zeros(numel(lambda),m);
Dispersion_material=zeros(numel(lambda),m);
M=zeros(N,N);
index=zeros(1,numel(rr));
actual_index=zeros(1,numel(rr));
index_profiles=zeros(numel(lambda),N); % records all the index profiles as the loops progresses
V_number=zeros(1,numel(lambda));
% opts for eigenvalue solver
clear opts

tic
ii=1;
for wavelength=[lambdastart:lambdastep:lambdaend]

% Calculate index of SiO2 at the current wavelength using eqn(6.11) 
% n_SiO2(ii)=sqrt(b1*(wavelength*10^6)^2/((wavelength*10^6)^2-a1)+b2*(wavelength*10^6)^2/((wavelength*10^6)^2-a2)+b3*(wavelength*10^6)^2/((wavelength*10^6)^2-a3)+1);
% Assume dispersion-less silica
    n_SiO2(ii)=1.4349;
% Construct Profile
    n_r1= sqrt(NA^2 + n_SiO2(ii)^2);
    n_r2= n_SiO2(ii);
 
%dont change ... building profile
    n_rr1=n_r1.*ones(1,numel(rr1));
    n_rr2=n_r2.*ones(1,numel(rr2));
    n_rr=cat(2,n_rr1,n_rr2);

% Build entire fiber profile for the specific wavelength
    actual_index = n_rr;
    
% Save profiles
    xls=[rr' actual_index'];

% Boundary conditions for matlab to do this computation - do not change  
    rrr=0; l=mod(L,2)+1;
    
    k=2*pi/wavelength;

%If the first radial position in the uploaded index profile file is zero,
%Matlab will divide by zero = bad situation! This is taken into account
%below.
xls(1,1) = del_r/10;
% Using the file loaded above, make an equally spaced index grid.
r = linspace(min(xls(:,1)),max(xls(:,1)),N);
index = interp1(xls(:,1),xls(:,2),r);
index_profiles(ii,:)=index;
%
% Build Matrix. This is the main part of the mode solver engine.
%
dr = r(2)-r(1);
%Inner part of diagonal matrix
for i = 2:N-1
	M(i,i-1) = 1/(dr^2)-1/(r(i)*2*dr);
	M(i,i)	= (index(i).^2)*k^2 - 2/dr.^2-(L^2)/(r(i)^2);
	M(i,i+1) = 1/(dr^2)+1/(r(i)*2*dr);
end 
%Outer part (boundaries) left part
	if l == 0 %field=0 @ r=0
	M(1,1) = (index(1).^2)*k^2 - 2/dr.^2-(L^2)/(r(1)^2);
	M(1,2) = 1/(dr^2)+1/(r(1)*2*dr);
	
    elseif l == 1 %field is symmetric @ r=0
	M(1,1) = (index(1).^2)*k^2 - 1/dr.^2-1/(r(1)*2*dr)-(L^2)/(r(1)^2);
	M(1,2) = 1/(dr^2)+1/(r(1)*2*dr);

    else %field is antisymmetric @ r=0
	M(1,1) = (index(1).^2)*k^2 - 3/dr.^2+1/(r(1)*2*dr)-(L^2)/(r(1)^2);
	M(1,2) = 1/(dr^2)+1/(r(1)*2*dr);
    end
    
%Right part
    if rrr == 0  %field=0 @ r: boundary
	M(N,N-1) = 1/(dr^2)-1/(r(N)*2*dr);
    M(N,N) = (index(N).^2)*k^2 - 2/dr.^2-(L^2)/(r(N)^2);
	
    elseif rrr == 1 %field is symmetric @ r: boundary
	M(N,N-1) = 1/(dr^2)-1/(r(N)*2*dr);
	M(N,N) = (index(N).^2)*k^2 - 1/dr.^2+1/(r(N)*2*dr)-(L^2)/(r(N)^2);

    else %field is antisymmetric @ r: boundary
	M(N,N-1) = 1/(dr^2)-1/(r(N)*2*dr);
	M(N,N) = (index(N).^2)*k^2 - 3/dr.^2-1/(r(N)*2*dr)-(L^2)/(r(N)^2);
    end


% Eigenvalues 'D' and Eigenvectors 'V' of diagonal matrix are found. 
M_sparse=sparse(M);
[V,D] = eigs(double(M_sparse),m,'lr');

% But the solutions found are all jumbled up. So, Sort the solutions.
% Largest beta (eigenvalue) first:
[beta_sqr,I] = sort(sum(D),'descend');
V = V(:,I);

% save variables
beta(ii,:)= sqrt(beta_sqr);
neff(ii,:)= sqrt(beta_sqr)*wavelength/2/pi;

%save the fields
for iii=1:1:m
   if (V(1,iii)<0)
        field(:,iii,ii)= -1* V(:,iii);
   elseif (V(1,iii)>0)
       field(:,iii,ii)= V(:,iii);
   end
end

V_number(ii)=k*r1*NA;

ii=ii+1;
end %end of wavelength loop
toc 

% Calculates dispersion
for modes=1:1:m
    nfit=fit(lambda',neff(:,modes),'cubicspline');
    dn=differentiate(nfit,lambda);
    dnfit=fit(lambda',dn,'cubicspline');
    d2neff=differentiate(dnfit,lambda);
    D=-lambda'.*d2neff/(c);
    Dispersion(:,modes)=D;
end


% Save fiber parameters into a .mat file
save(['data.mat'],'beta','lambda','neff','field','r','index_profiles','N','r1','r2','L','m','V_number','Dispersion','w')


% Explanation of variables

% all in SI units

% beta(lambda_indices,mode): beta values (1/m)
% lambda: all the discrete wavelengths that waveguide is simulated for (m)
% neff(lambda_indices,mode): effective index of particular mode and wavelength
% field(:,mode,lambda_indices): fields of the mode vs r for a particular mode and lambda_indices
% r: radius in steps of (total_radius/N) (m)
% index_profiles(lambda_indices,:): records the actual index (including material dispersion) as a function of r
% N: N is the number of grid points the mesh is made over
% V_number(lambda_indices): V number of a fiber
% Dispersion(lambda_indices,mode): dispersion (s/m.m)

% END OF PROGRAM

%% PLOTTING
%====================================================================
% Section A
%
% The eignvectors are given in the tensor field. The first index 
% corresponds to the radial distribution of the field, for a given mode
% order m (second index of tensor), and at a given wavelength, as specified
% in the wavelength loop (third index of the tensor). Note that the
% wavelength index of this tensor is the same as ii, in the the for-do loop
% above. Plotting the first index (axis) of the tensor versus radial
% position r gives the *radial* field distribution.
%
%Below is an illustrative example. We ask matlab to compute beta and neff
%for the mth mode with angular index L , and print it along with the plot 
% of the radial distribution of this mode. (This plot also shows the 
% refractive index profile to identify the position of the mode with respect 
%to the waveguide.

load('data.mat')

m=1;

input_wavelength=1.0e-6;
input_index=round((input_wavelength-lambda(1))/(lambda(2)-lambda(1)))+1;


figure(1); 
[AX,H1,H2] = plotyy(r*10^6,index_profiles(input_index,:), r*10^6,field(:,m,input_index).^2);
set(get(AX(1),'Ylabel'),'String','Refractive Index')
set(get(AX(2),'Ylabel'),'String','Mode Intensity (a.u.)')
xlabel('Radial Position (um)'); 

title(['LP(' num2str(L) ',' num2str(m), ') mode at '...
    num2str(lambda(input_index)*10^9) 'nm;  \beta^2='...
    num2str((beta(input_index,m)*10^-6)^2) '/um^2' '  =>  \beta='...
    num2str(beta(input_index,m)*10^-6) '/um' '  =>   neff='...
    num2str(neff(input_index,m))])

%====================================================================
%SECTION B
%
% The plot above only visualises the radial field. But if there is an
% angular component to it too (for e.g., for all modes with L>0), the polar
% plot, given below, provides a better visualisation. Note that the fields
% are squared - so, this plot is that of mode intensity, not mode field.
%
theta = linspace(0,2*pi,N);
[th,R] = meshgrid(theta,r);
[X,Y] = pol2cart(th,R);
figure(2)
mesh(X,Y,(field(:,m,input_index)*cos(L*theta)).^2);
view(45,75)

%========================================================================
% Section C
% In this section, we plot the wavelength dependent neff and dispersion. 
% Note that mode order m, defined in section A, is used here too.

figure(3)
plot(lambda*10^9,neff(:,:))
hold on; plot(lambda*10^9, index_profiles(:,N), 'linewidth',2)
xlabel('Wavelength (nm)')
ylabel('Effective Index')

figure(4)
plot(lambda*10^9,(10^-9*10^3)/10^-12*(Dispersion(:,:)))
xlabel('Wavelength (nm)')
ylabel('Dispersion (ps/nm-km)')