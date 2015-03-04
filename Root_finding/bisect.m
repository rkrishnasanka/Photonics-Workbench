function [root,yval] = bisect(varfun, a, b, tol)

%We first check to see if there is the needed sign change at 
%the two function values
ya = feval(varfun, a);
yb = feval(varfun, b);

if sign(ya) == sign(yb)
    error('Function has same sign at endpoints')
end

%We assign the default tolerance, if none is specified
if nargin <4
    tol = eps*max([abs(a) abs(b) 1]);
end

%We initialize the iteration

an = a; %Limit a after n iterations
bn = b; %Limit b after n iterations
n = 0; %Iteration count

%Bisection Loop
while (b-a)/2^n >= tol
    xn = (an +bn)/2; %Find midpoint of the limits
    yn = feval(varfun, xn); %Find the solution of the function at the midpoint
    n = n+1; %Increase the iteration by 1 
    
    if yn == 0
        fprintf('numerically exact root')
        
        root = xn;
        yval = yn;
        return
    elseif sign(yn) == sign(ya)
        an = xn;
        ya = yn;
    else
        bn = xn;
        yb = yn;
    end
end

root = xn;
yval = yn;
    