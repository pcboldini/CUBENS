 
function [C,kappa,Pr,Rho,T] = Calc_PR_Prop(h_ref,PR_PROP,BL_PARA)

global visc

Zc=PR_PROP.Zc;
Zc_1=PR_PROP.Zc^(-1);
Zc_2=PR_PROP.Zc^(-2);
a=PR_PROP.a;
b=PR_PROP.b;
Ru=PR_PROP.Ru;
K=PR_PROP.K;

CvR=BL_PARA.CvR;
p_inf=BL_PARA.p_inf;
Rho_inf=BL_PARA.Rho_inf;
Cp_inf=BL_PARA.Cp_inf;
Pr_inf=BL_PARA.Pr_inf;
h_inf=BL_PARA.h_inf;
T_inf=BL_PARA.T_inf;
mu_inf=BL_PARA.mu_inf;
kappa_inf=BL_PARA.kappa_inf;

T0=interp1(PR_PROP.h0_approx,PR_PROP.T0_approx,h_ref*h_inf,'spline');
rho0=interp1(PR_PROP.T0_approx,PR_PROP.rho0_approx,T0,'spline');
v0=1/rho0;

fun = @PR2d;
x0 = [T0,v0];
options = optimoptions('fsolve','Display','none');
x = fsolve(fun,x0,options);
T=x(1);
v=x(2);
Rho=1/v;
alpha=(1+K*(1-sqrt(T)))^2;

if strcmp(visc,'VDW')
    
    if T<=1.50
        mu_1=34*10^(-5)*T^(0.94);
    else
        mu_1=17.78*10^(-5)*(4.58*T-1.67)^(5/8);
    end 
 
   f_rho=0.10230+0.023364*Rho+0.058533*Rho.^2-0.040758*Rho.^3+0.0093324*Rho.^4;
   mu_diff=(f_rho.^4-10^(-4));
   mu=mu_diff+mu_1;

   mu=mu/mu_inf;
   
   kappa_RT_tri=6.54*10^(-5).*(exp(0.2826.*T)-1./exp(0.3976.*T.^2));
   
   kappa_factor=(0.307*CvR+0.539);

   kappa_EU=15/4*Ru*mu_1*kappa_factor;
   
   if Rho<0.50
      f_rho=14.0*(exp(0.535*Rho)-1);
      kappa_diff=(f_rho*10^(-8))/Zc^5;
      kappa=(kappa_diff)*(4.1868/10^(-2))+kappa_EU;
   elseif Rho>=0.50 && Rho<2.0
       f_rho_2=13.1*(exp(0.67*Rho)-1.069);
       kappa_diff=(f_rho_2*10^(-8))/Zc^5;
       kappa=(kappa_diff)*(4.1868/10^(-2))+kappa_EU;
   else
       f_rho_3=2.976*(exp(1.155*Rho)+2.016); 
       kappa_diff=(f_rho_3*10^(-8))/Zc^5;
       kappa=(kappa_diff)*(4.1868/10^(-2))+kappa_EU;
   end
         
   kappa=kappa/kappa_inf;

elseif strcmp(visc,'Chung') 
    [mu,kappa] = Chung(Rho,T,PR_PROP);
    mu=mu/mu_inf;
    kappa=kappa/kappa_inf;   
   
elseif strcmp(visc,'constant')    
    
    mu=1;
    kappa   = mu; % Thermal conductivity

end

Cv_rho=CvR-(a*K*(K+1))/(4*b*sqrt(2*T))*log((1+b*Zc_1*Rho*(1-sqrt(2)))/(1+b*Zc_1*Rho*(1+sqrt(2))));
dPdT_Rho=K*sqrt(alpha/T)*(Rho^2*a*Zc_2)/(1+2*b*Zc_1*Rho-b^2*Zc_2*Rho^2)-(Rho*Zc_1)/(Rho*b*Zc_1-1);
dPdRho_T=Zc_1*T/(Zc_1*b*Rho-1)^2-a*Zc_2*alpha*(2*Rho+2*b*Zc_1*Rho^2)/(1+2*b*Zc_1*Rho-Rho^2*b^2*Zc_2)^2;
Cp_rho=Cv_rho+T/Rho^2*Zc*dPdT_Rho^2/dPdRho_T;

Pr=Cp_rho*mu/kappa*Pr_inf/Cp_inf;
C=Rho/Rho_inf*mu;

Rho=Rho/Rho_inf;
T=T/T_inf;

function F = PR2d(x) %%% x(1)=T, x(2)=v

F(1) = x(1)*Zc_1/(x(2)-b*Zc_1) - a*((1+K*(1-sqrt(x(1))))^2)*Zc_2/(x(2)*(x(2)+2*b*Zc_1)-b^2*Zc_2) - p_inf;
F(2) = CvR*x(1)*Zc_1+a*Zc_1/(2*sqrt(2)*b)*((K+1)^2-K*(K+1)*sqrt(x(1)))*log((1+b*Zc_1/x(2)*(1-sqrt(2)))/(1+b*Zc_1/x(2)*(1+sqrt(2))))+p_inf*x(2)-h_ref*h_inf;

end

end
