function n_mesh=refindex(x,interface,index_layer)
%assign the values of refractive indices to mesh points
N_mesh=length(x);
NumberOfLayers=length(index_layer);
i_mesh=1;
for k=1:NumberOfLayers      %loop over all layers
    for i=1:interface(k)    %loop within layer
        n_mesh(i_mesh+1)=index_layer(k);
        i_mesh=i_mesh+1;
    end
end
