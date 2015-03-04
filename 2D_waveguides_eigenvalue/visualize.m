function visualize(r,m,fn)
close;
[~,h0]=contour(r.y,r.x,r.mode(:,:,m),48);
hold on;
h1=plot([r.CW(3),r.CW(4)],[r.RB(1),r.RB(1)]);
h2=plot([r.RB(3),r.RB(3)],[r.RB(1),r.RB(2)]);
h3=plot([r.RB(4),r.RB(4)],[r.RB(1),r.RB(2)]);
h4=plot([r.RB(3),r.RB(4)],[r.RB(2),r.RB(2)]);
hold off;
axis equal;
set([h0,h1,h2,h3,h4],'LineWidth',1.5);
if ~strcmp(fn,'')
epsfn=[fn,'.eps'];
print('-depsc',epsfn);
command=['bash ','epstopdf ',epsfn];
system(command);
delete(epsfn);
end;
end % function visualize