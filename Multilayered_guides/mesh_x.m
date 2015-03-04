function x=mesh_x(d_s,d_layer,d_c,NumberMesh)
%generate 1D mesh along x

d_total=[d_s,d_layer,d_c];          %thickness of all layers
NumberOfLayers=length(d_total);     %determine of all layers
delta=d_total./NumberMesh;
x(1)=0.0;                           %coordinate of first mesh point
i_mesh=1;
for k=1:NumberOfLayers             %loop over all layers
    for i=1:NumberMesh(k)          %loop within layer
        x(i_mesh+1)=x(i_mesh)+delta(k);
        i_mesh=i_mesh+1;
    end
end

   