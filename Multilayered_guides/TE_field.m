function TE_mode_field=TE_field(beta,index_mesh,x,k_zero)
%determine TE optical field for all layers
TotalMesh=length(x);
zz=beta*k_zero;
kappa=0;
for n=1:(TotalMesh)
    kappa(n)=sqrt((k_zero*index_mesh(n))^2-zz^2);
end

U(1)=1.0;
temp=imag(kappa(1));
if(temp<0), kappa(1)=-kappa(1);
end
V(1)=kappa(1);

for n=2:(TotalMesh)
    cc=cos(kappa(n)*(x(n)-x(n-1)));
    ss=sin(kappa(n)*(x(n)-x(n-1)));
    m(1,1)=cc;
    m(1,2)=-1i/kappa(n)*ss;
    m(2,1)=-1i*kappa(n)*ss;
    m(2,2)=cc;
    
    U(n)=m(1,1)*U(n-1)+m(1,2)*V(n-1);
    V(n)=m(2,1)*U(n-1)+m(2,2)*V(n-1);
    
end
TE_mode_field=abs(U);
max_value=max(TE_mode_field);
h=plot(x,TE_mode_field/max_value);
xlabel('x (microns)','FontSize',22);
ylabel('TE electric field','FontSize', 22);
set(h,'LineWidth',1.5);
set(gca,'FontSize',22);
pause
close all
    
    
    