function rwg=solve(r)
% calculate 6 eigenvectors of sparse matrix
% algebraically largest eigenvalues
[efun,eval]=eigs(r.hh,6,'la');
eval=diag(eval);
% guided modes
guided=(r.ES<eval);
modes=efun(:,guided);
[NX,NY]=size(r.eps);
rwg=r;
% modes must be reshaped
rwg.mode=reshape(modes,[NY,NX,sum(guided)]);
end % function solve