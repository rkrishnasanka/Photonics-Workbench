%Driver file to calculate propagation constants and electric field
%profile (TE) of arbitrary multilayered waveguides

clear all;
format long;

%input structure for analysis 
structure;
%
epsilon=1e-6;
TE_mode=[];
n_max=max(n_layer);
z1=n_max; %max value of refr index
n_min=max(n_s,n_c)+0.001; %min value of refr index
dz=0.005; %iteration step
mode_control=0;
%
while(z1>n_min)
    
    z0=z1-dz; %startin point for Muller method
    z2=0.5*(z1+z0);
    z_new = muller(@f_TE , z0, z1, z2);
    if (z_new ~= 0)         %verifying mode existence
        
        for u=1:length(TE_mode)
            if(abs(TE_mode(u)-z_new) < epsilon)
                mode_control=1; break; %mode found
            end
        end
        
        if (mode_control ==1)
            mode_control=0;
        else
            TE_mode(length(TE_mode)+1)=z_new;
        end
    end
    z1=z0;
end
%

TE_mode=sort(TE_mode, 'descend');
TE_mode'    %outputs all calculated modes
beta=TE_mode(1); %select mode for plotting field profile
x=mesh_x(d_s,d_layer,d_c,NumberMesh);
n_total=[n_s, n_layer,n_c];
n_mesh=refindex(x,NumberMesh,n_total);
TE_mode_field=TE_field(beta,n_mesh,x,k_0);

