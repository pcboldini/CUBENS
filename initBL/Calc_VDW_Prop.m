 
function [C,kappa,Pr,Rho,T] = Calc_VDW_Prop(h_ref,VDW_PROP,BL_PARA)

global visc

Zc=VDW_PROP.Zc;
R=VDW_PROP.R;
a=VDW_PROP.a;
Ru=VDW_PROP.Ru;

CvR=BL_PARA.CvR;
p_inf=BL_PARA.p_inf;
Rho_inf=BL_PARA.Rho_inf;
Cp_inf=BL_PARA.Cp_inf;
Pr_inf=BL_PARA.Pr_inf;
h_inf=BL_PARA.h_inf;
T_inf=BL_PARA.T_inf;
mu_inf=BL_PARA.mu_inf;
kappa_inf=BL_PARA.kappa_inf;

T0=interp1(VDW_PROP.h0_approx,VDW_PROP.T0_approx,h_ref*h_inf,'spline');
rho0=interp1(VDW_PROP.T0_approx,VDW_PROP.rho0_approx,T0,'spline');

fun = @VDW2d;
x0 = [T0,rho0];
options = optimoptions('fsolve','Display','none');
x = fsolve(fun,x0,options);
T=x(1);
Rho=x(2);

if strcmp(visc,'JST')
   
   [mu,kappa] = JST(Rho,T,BL_PARA,VDW_PROP); %% non-dimensional
   
   mu=mu/mu_inf;
   kappa=kappa/kappa_inf;
   
elseif strcmp(visc,'Constant')    
    
    mu=1;
    kappa   = mu; % Thermal conductivity

elseif strcmp(visc,'Chung') 

    [mu,kappa] = Chung(Rho,T,VDW_PROP);
    mu=mu/mu_inf;
    kappa=kappa/kappa_inf;

end

sol_v=1/Rho;
Cp=CvR+1/(1-(3*sol_v-1).^2/(4*T*sol_v.^3));

Pr=Cp*mu/kappa*Pr_inf/Cp_inf;
C=Rho/Rho_inf*mu;

Rho=Rho/Rho_inf;
T=T/T_inf;

function F = VDW2d(x) %%% x(1)=T, x(2)=rho
    
F(1) = (8*x(1))/(3*(1/x(2))-1) - 3*x(2)^2 - p_inf ;
F(2) = R*CvR*x(1)-a*x(2) +p_inf/x(2) - h_ref*h_inf;

end

end
