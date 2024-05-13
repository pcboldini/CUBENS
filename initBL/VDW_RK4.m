% blasius compressible similarity solution

function [fSol,CSol,KappaSol,PrSol,RhoSol,TSol]=VDW_RK4(BL_PARA,VDW_PROP,N)

global wall_bc

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

T_inf=BL_PARA.T_inf;

% Wall Temperature: Default is the recovery temperature:
T_rec = ( 1+sqrt(BL_PARA.Pr_inf)*(BL_PARA.gamma-1)/2*BL_PARA.M_inf*BL_PARA.M_inf )* T_inf;
[h_rec,~,~,~,~,~,~,~] = VDW(T_rec,VDW_PROP,BL_PARA);
h_inf=BL_PARA.h_inf;

%% Shooting solver for BL solution

% Tolerance of Newton-Raphson Method 
tol=1e-13;

% Numerical delta for Jacobian derivatives
epsilon=1e-10;

if strcmp(wall_bc,'adiab')
    
%    [ Cf''   g ]
ic = [0.5    h_rec/h_inf ];

% Shooting method for initial condition
fprintf('Iteration process: adiabatic wall (Ny=%i) \n',N)
[eta,fSol,CSol,KappaSol,PrSol,RhoSol,TSol,Initialcond] = fsCompAd(ic,BL_PARA,VDW_PROP,N,dh);

elseif strcmp(wall_bc,'isoth')
    
% Enthalpy at the wall

[h_wall,~,~,~,~,~,~,~] = VDW(BL_PARA.T_wall,VDW_PROP,BL_PARA);
    
%    [ Cf''   g' ]
ic = [0.5    0 ];

g_w=h_wall/h_inf;

% Shooting method for initial condition
[eta,fSol,CSol,KappaSol,PrSol,RhoSol,TSol,Initialcond] = fsCompIso(ic,g_w,BL_PARA,VDW_PROP,N,dh);

end

function [eta,fSol,CSol,KappaSol,PrSol,RhoSol,TSol,v] = fsCompAd(v,BL_PARA,VDW_PROP,N,dh)

%FSCOMP compressible boundary layer
% computation of similarity solution
% by shooting method

% input values
Pr_inf=BL_PARA.Pr_inf;
M_inf=BL_PARA.M_inf;
Ec_inf=BL_PARA.Ec_inf;
Cp_inf=BL_PARA.Cp_inf;
h_inf=BL_PARA.h_inf;
T_inf=BL_PARA.T_inf;
Zc=VDW_PROP.Zc;
ue2he =Ec_inf*Cp_inf*T_inf/h_inf/Zc;
countit=1;
count=0;
    
v = newtonRaphson2(@residual,v,tol,epsilon);

    function r = residual(v)
        r = zeros(length(v),1);
        [xSol,ySol,CSol,KappaSol,PrSol,RhoSol,TSol] = VDW_RK4_BF(@fsCompOde,inCond(v),VDW_PROP,BL_PARA,N,dh);
        % Check
        count=count+1;
        lastRow = size(ySol,1);
        r(1) = ySol(lastRow,2) - 1;  % f'-1
        r(2) = ySol(lastRow,4) - 1;  % g -h_inf
        if numel(v)+1==count/countit
        fprintf ('Iteration %i: T_w = %.6f, Cf''(end) = %.6f, g(end) = %.6f \n',countit,TSol(1),ySol(lastRow,2),ySol(lastRow,4));
        countit=countit+1;
        end
    end
    %----------------------------------------------------------------------
    
        % nested functions
    function dfdeta = fsCompOde(eta,f,C,Pr,Rho)    
        dfdeta = zeros(1,6);
        dfdeta(1) =  f(2);
        dfdeta(2) =  f(3)/C;
        dfdeta(3) = -f(3)*f(1)/C;
        dfdeta(4) =  Pr*f(5)/C;
        dfdeta(5) = -Pr*f(1)*f(5)/C - ue2he*f(3)^2/C;
        dfdeta(6) =  sqrt(2)/Rho;
    end
    % ---------------------------------------------------------------------

    function f = inCond(V)
        %   [f f' Cf''  g   Cg' eta];
        f = [0  0  V(1)  V(2) 0   0];
    end
    %----------------------------------------------------------------------
    
fprintf ('Iteration completed! Inital solution at the wall is: Cf''=%.4f, g=%.4f \n',v(1),v(2));

% Calculation of the BL with the correct BC at the wall
[eta,fSol,CSol,KappaSol,PrSol,RhoSol,TSol] = VDW_RK4_BF(@fsCompOde,inCond(v),VDW_PROP,BL_PARA,N,dh);


end % fsComp

function [eta,fSol,CSol,KappaSol,PrSol,RhoSol,TSol,v] = fsCompIso(v,g_w,BL_PARA,VDW_PROP,N,dh)

%FSCOMP compressible boundary layer
% computation of similarity solution
% by shooting method

% input values
Pr_inf=BL_PARA.Pr_inf;
M_inf=BL_PARA.M_inf;
Ec_inf=BL_PARA.Ec_inf;
Cp_inf=BL_PARA.Cp_inf;
h_inf=BL_PARA.h_inf;
T_inf=BL_PARA.T_inf;
Zc=VDW_PROP.Zc;
ue2he =Ec_inf*Cp_inf*T_inf/h_inf/Zc;
countit=1;
count=0;
    
%% Iteration Method
% Newton-Raphson 2D 
v = newtonRaphson2(@residual,v,tol,epsilon);

    function r = residual(v)
        r = zeros(length(v),1);
        [xSol,ySol,CSol,KappaSol,PrSol,RhoSol,TSol] = VDW_RK4_BF(@fsCompOde,inCond(v),VDW_PROP,BL_PARA,N,dh);
        % Check
        count=count+1;
        lastRow = size(ySol,1);
        r(1) = ySol(lastRow,2) - 1;  % f'-1
        r(2) = ySol(lastRow,4) - 1;  % g -1
        if numel(v)+1==count/countit
        fprintf ('Iteration %i: dhdy_w = %.6f, Cf''(end) = %.6f, g(end) = %.6f \n',countit,ySol(1,5),ySol(lastRow,2),ySol(lastRow,4));
        countit=countit+1;
        end
    end
    %----------------------------------------------------------------------
    
        % nested functions
    function dfdeta = fsCompOde(eta,f,C,Pr,Rho)    
        dfdeta = zeros(1,6);
        dfdeta(1) =  f(2);
        dfdeta(2) =  f(3)/C;
        dfdeta(3) = -f(3)*f(1)/C;
        dfdeta(4) =  Pr*f(5)/C;
        dfdeta(5) = -Pr*f(1)*f(5)/C - ue2he*f(3)^2/C;
        dfdeta(6) =  sqrt(2)/Rho;
    end
    % ---------------------------------------------------------------------

    function f = inCond(V)
        %   [f f' Cf''  g   Cg' eta];
        f = [0  0  V(1) g_w  V(2)   0];
    end
    
        fprintf ('Iteration completed! Inital solution at the wall is: Cf''=%.4f, g''=%.4f \n',v(1),v(2));

    %----------------------------------------------------------------------
% Calculation of the BL with the correct BC at the wall
[eta,fSol,CSol,KappaSol,PrSol,RhoSol,TSol] = VDW_RK4_BF(@fsCompOde,inCond(v),VDW_PROP,BL_PARA,N,dh);

end % fsComp

 
end