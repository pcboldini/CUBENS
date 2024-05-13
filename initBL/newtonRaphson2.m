function root = newtonRaphson2(func,x,tol,epsilon)
% Newton-Raphson method of finding a root of simultaneous
% equations fi(x1,x2,...,xn) = 0, i = 1,2,...,n.
% USAGE: root = newtonRaphson2(func,x,tol)
% INPUT:
% func = handle of function that returns[f1,f2,...,fn].
% x = starting solution vector [x1,x2,...,xn].
% tol = error tolerance (default is 1.0e4*eps).
% OUTPUT:
% root = solution vector.

if size(x,1) == 1; x = x'; end
% x must be column vector
for i = 1:15
    [jac,f0] = Jacobian(func,x,epsilon);
    if sqrt(dot(f0,f0)/length(x)) < tol
        root = x; return
    end
    dx = jac\(-f0);
    x = x + dx;
    error=sqrt(dot(dx,dx)/length(x));
    fprintf ('Error: %d \n',error);
    if error < tol*max(abs(x),1.0)
        root = x; return
    end
end
root = x; return
end

