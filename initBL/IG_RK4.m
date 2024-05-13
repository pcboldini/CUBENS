% blasius compressible similarity solution

function [fSol]=IG_RK4(BL_PARA,N)

global wall_bc visc

%% Mesh
ymax = 10;
yi=3;

% Chebyshev Configuration Point
z = -cos((0:(N-1))/(N-1)*pi);

a = ymax * yi / (ymax - 2.0 * yi);
b = 1 + 2.0 * a / ymax;  

    for ii=1:N
        y(ii) = a * (1 + z(ii)) / (b - z(ii));
    end

for kk=1:numel(y)-1
dh(kk)=y(kk+1)-y(kk);
end


%% Shooting solver for BL solution

% Tolerance of Newton-Raphson Method 
tol=1e-13;

% Numerical delta for Jacobian derivatives
epsilon=1e-10;

% Wall Temperature: Default is the recovery temperature:
T_rec = ( 1+sqrt(BL_PARA.Pr_inf)*(BL_PARA.gamma-1)/2*BL_PARA.M_inf*BL_PARA.M_inf ) ;

if strcmp(wall_bc,'adiab')
    
%    [ Cf''   g ]
ic = [.05    T_rec ];

% Shooting method for initial condition
[eta,fSol,Initialcond] = fsCompAd(ic,BL_PARA,N,dh);

elseif strcmp(wall_bc,'isoth')
    
% Approximation of the enthalpy gradient at the wall

T_star=0.5*BL_PARA.T_wall+0.22*T_rec+0.28*BL_PARA.T_inf;
F_c=T_star/BL_PARA.T_inf;
rho_w = BL_PARA.p_inf/(BL_PARA.R_inf*BL_PARA.T_wall);
[~,kappa_w]=Sutherland_dim(BL_PARA);

St=BL_PARA.Pr_inf^(-2/3)*0.664/(2*F_c);

dhdy_w=sqrt(2)*St*BL_PARA.rho_inf/rho_w*BL_PARA.mu_inf*BL_PARA.cp_inf^2/kappa_w*(T_rec-BL_PARA.T_wall)/BL_PARA.h_inf;

%    [ Cf''   g' ]
ic = [.05    0 ];
g_w=BL_PARA.T_wall/BL_PARA.T_inf;

% Shooting method for initial condition
[eta,fSol,Initialcond] = fsCompIso(ic,g_w,BL_PARA,N,dh);

end


function [eta,fSol,v] = fsCompAd(v,BL_PARA,N,dh)

%FSCOMP compressible boundary layer
% computation of similarity solution
% by shooting method

% parameters
gamma = BL_PARA.gamma;
S    = BL_PARA.S_mu_ref;
powerexp    = BL_PARA.powerexp;
M_inf    = BL_PARA.M_inf;
T_inf = BL_PARA.T_inf;

% input values
ue2he = (gamma-1)*M_inf^2;
TsTe  = S/T_inf;
Pr_inf    = BL_PARA.Pr_inf;
countit=1;
count=0;
    
v = newtonRaphson2(@residual,v,tol,epsilon);

[eta,fSol] = RK4(@fsCompOde,inCond(v),N,dh);

    function r = residual(v)
        %global XSTART XSTOP H
        r = zeros(length(v),1);
        [xSol,ySol] = RK4(@fsCompOde,inCond(v),N,dh);
        % Check
        count=count+1;
        lastRow = size(ySol,1);
        r(1) = ySol(lastRow,2) - 1;  % f'-1
        r(2) = ySol(lastRow,4) - 1;  % g -1
        if numel(v)+1==count/countit
%        fprintf ('Iteration %i: T_w = %.6f, Cf''(end) = %.6f, g(end) = %.6f \n',countit,ySol(1,4),ySol(lastRow,2),ySol(lastRow,4));
        countit=countit+1;
        end
    end
    %----------------------------------------------------------------------
    
        % nested functions
    function dfdeta = fsCompOde(eta,f)    
        dfdeta = zeros(1,6);
         if strcmp(visc,'Constant')
             C = 1./f(4);
         elseif  strcmp(visc,'PowerLaw')
             C = f(4)^(powerexp-1);
         elseif  strcmp(visc,'Sutherland')
             C = sqrt(f(4)).*(1+TsTe)./(f(4)+TsTe);
         end
        dfdeta(1) =  f(2);
        dfdeta(2) =  f(3)/C;
        dfdeta(3) = -f(3)*f(1)/C;
        dfdeta(4) =  Pr_inf*f(5)/C;
        dfdeta(5) = -Pr_inf*f(1)*f(5)/C - ue2he*f(3)^2/C;
        dfdeta(6) =  sqrt(2)*f(4);
    end
    % ---------------------------------------------------------------------

    function f = inCond(V)
        %   [f f' Cf''  g   Cg' eta];
        f = [0  0  V(1)  V(2) 0   0];
    end
    %----------------------------------------------------------------------
end % fsComp

function [eta,fSol,v] = fsCompIso(v,g_w,BL_PARA,N,dh)

%FSCOMP compressible boundary layer
% computation of similarity solution
% by shooting method

% parameters
gamma = BL_PARA.gamma;
S    = BL_PARA.S_mu_ref;
powerexp    = BL_PARA.powerexp;
M_inf    = BL_PARA.M_inf;
T_inf = BL_PARA.T_inf;

% input values
ue2he = (gamma-1)*M_inf^2;
TsTe  = S/T_inf;
Pr_inf    = BL_PARA.Pr_inf;
countit=1;
count=0;
    
v = newtonRaphson2(@residual,v,tol,epsilon);

[eta,fSol] = RK4(@fsCompOde,inCond(v),N,dh);

    function r = residual(v)
        %global XSTART XSTOP H
        r = zeros(length(v),1);
        [xSol,ySol] = RK4(@fsCompOde,inCond(v),N,dh);
        % Check
        count=count+1;
        lastRow = size(ySol,1);
        r(1) = ySol(lastRow,2) - 1;  % f'-1
        r(2) = ySol(lastRow,4) - 1;  % g -1
        if numel(v)+1==count/countit
 %       fprintf ('Iteration %i: dhdy_w = %.6f, Cf''(end) = %.6f, g(end) = %.6f \n',countit,ySol(1,5),ySol(lastRow,2),ySol(lastRow,4));
        countit=countit+1;
        end
    end
    %----------------------------------------------------------------------
    
        % nested functions
    function dfdeta = fsCompOde(eta,f)    
        dfdeta = zeros(1,6);
         if strcmp(visc,'Constant')
             C = 1./f(4);
         elseif  strcmp(visc,'PowerLaw')
             C = f(4)^(powerexp-1);
         elseif  strcmp(visc,'Sutherland')
             C = sqrt(f(4)).*(1+TsTe)./(f(4)+TsTe);
         end
        dfdeta(1) =  f(2);
        dfdeta(2) =  f(3)/C;
        dfdeta(3) = -f(3)*f(1)/C;
        dfdeta(4) =  Pr_inf*f(5)/C;
        dfdeta(5) = -Pr_inf*f(1)*f(5)/C - ue2he*f(3)^2/C;
        dfdeta(6) =  sqrt(2)*f(4);
    end
    % ---------------------------------------------------------------------

    function f = inCond(V)
        %   [f f' Cf''  g   Cg' eta];
        f = [0  0  V(1) g_w  V(2)   0];
    end
    %----------------------------------------------------------------------
end % fsComp


    
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
for i = 1:20
    [jac,f0] = jacobian(func,x,epsilon);
    if sqrt(dot(f0,f0)/length(x)) < tol
        root = x; return
    end
    dx = jac\(-f0);
    x = x + dx;
    error=sqrt(dot(dx,dx)/length(x));
%    fprintf ('Error: %d \n',error);
    if error < tol*max(abs(x),1.0)
        root = x; return
    end
end
root = x; return
end

function [jac,f0] = jacobian(func,x,h)
% Returns the Jacobian matrix and f(x).

n = length(x);
jac = zeros(n);
f0 = feval(func,x);
for i =1:n
    temp = x(i);
    x(i) = temp + h;
    f1 = feval(func,x);
    x(i) = temp;
    jac(:,i) = (f1 - f0)/h;
end
end

function[xSol,ySol] = RK4(dEqs,y,N,dh)
x=0;
xSol = zeros(N,numel(x));
ySol = zeros(N,numel(y));
xSol(1) = x; ySol(1,:) = y;
k = 1;
for i=1:N-1
        k1=feval(dEqs,x,y);
        k2=feval(dEqs,x+dh(i)/2,y+k1*dh(i)/2);
        k3=feval(dEqs,x+dh(i)/2,y+k2*dh(i)/2);
        k4=feval(dEqs,x+dh(i),y+k3*dh(i));
        dy=(dh(i)/6)*(k1+2*k2+2*k3+k4);
        x = x+dh(i);
        y = y + dy;
        k = k+1;
        xSol(k) = x; ySol(k,:) = y;  
end
end

end