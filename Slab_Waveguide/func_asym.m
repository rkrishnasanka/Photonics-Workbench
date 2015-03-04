function f=func_asym(beta,n_c,n_s,n_f,k,h)

%constructs the search function for the asymmetric 3 layer slab

gamma_s=sqrt(beta.^2-(n_s*k)^2);
gamma_c=sqrt(beta.^2-(n_c*k)^2);
kappa_f=sqrt((n_f*k)^2-beta.^2);

denom=kappa_f-(gamma_c.*gamma_s)./kappa_f;
f=tan(kappa_f*h)-(gamma_s+gamma_c)./denom;