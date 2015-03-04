function rwg=digitize(r,NX,NY)
x=linspace(r.CW(1),r.CW(2),NX);
y=linspace(r.CW(3),r.CW(4),NY);
[X,Y]=meshgrid(x,y);
RIB=(X>r.RB(1))&(X<r.RB(2))&(Y>r.RB(3))&(Y<r.RB(4));
SUB=(X<=r.RB(1));
COV=(X>r.RB(1))&~RIB;
prm=r.EC*COV+r.ES*SUB+r.ER*RIB;
rwg=r;
rwg.x=x;
rwg.y=y;
rwg.eps=prm;
% note that NX and NY may be reconstructed
% from rwg.eps
end % function digitize