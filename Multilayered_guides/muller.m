%this function implements the Muller's root finding method for complex
%roots
function f_val=muller(f, x0, x1, x2)

iter_max=100;
f_tol=1e-6;
x_tol=1e-6;
y0=f(x0);
y1=f(x1);
y2=f(x2);
iter=0;
while(iter<=iter_max)
    iter=iter+1;
    a=( (x1-x2)*(y0-y2) - (x0-x2)*(y1-y2))/( (x0-x2)*(x1-x2)*(x0-x1) );
    %
    b=( (x0-x2)^2*(y1-y2) - (x1-x2)^2*(y0-y2))/( (x0-x2)*(x1-x2)*(x0-x1) );
    %
    c=y2;
    %
    if (a~=0)
        D=sqrt(b*b-4*a*c);
        q1=b+D;
        q2=b-D;
        if (abs(q1)<abs(q2))
            dx=-2*c/q2;
        else
            dx=-2*c/q1;
        end
    elseif (b~=0)
        dx=-c/b;
    else
        warning('Muller method failed to find a root')
        break;
    end
    x3=x2 +dx;
    x0=x1;
    x1=x2;
    x2=x3;
    y0=y1;
    y1=y2;
    y2=feval(f,x2);
    if (abs(dx)<x_tol && abs(y2)<f_tol)
        break;
    end
end
%check that only proper values are calcualated
if (abs(y2)<f_tol)
    f_val=x2;
    return;
else
    f_val=0;
end

