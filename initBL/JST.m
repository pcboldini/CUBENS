function [mu_VDW,kappa_VDW] = JST(rho,T,BL_PARA,VDW_PROP)
%JST Summary of this function goes here
%   Detailed explanation goes here

% Loading
CvR=BL_PARA.CvR;
Ru=VDW_PROP.Ru;
Zc=VDW_PROP.Zc;
kappa_factor=(0.307*CvR+0.539);

% Viscosity

if T<=1.50
    mu_1=34*10^(-5)*T^(0.94);
else
    mu_1=17.78*10^(-5)*(4.58*T-1.67)^(5/8);
end
 
f_rho=0.10230+0.023364*rho+0.058533*rho.^2-0.040758*rho.^3+0.0093324*rho.^4;
mu_diff=(f_rho.^4-10^(-4));
mu_VDW=mu_diff+mu_1;

% Thermal conductivity      
        
kappa_EU=15/4*Ru*mu_1*kappa_factor;

kappa_RT_tri=6.54*10^(-5).*(exp(0.2826.*T)-1./exp(0.3976.*T.^2));  
     
if rho<0.50
    f_rho=14.0*(exp(0.535*rho)-1);
    kappa_diff=(f_rho*10^(-8))/Zc^5;
    kappa_VDW=kappa_diff*(4.1868/10^(-2))+kappa_EU;
elseif rho>=0.50 && rho<2.0
    f_rho_2=13.1*(exp(0.67*rho)-1.069);
    kappa_diff=(f_rho_2*10^(-8))/Zc^5;
    kappa_VDW=kappa_diff*(4.1868/10^(-2))+kappa_EU;
else
    f_rho_3=2.976*(exp(1.155*rho)+2.016); 
    kappa_diff=(f_rho_3*10^(-8))/Zc^5;
    kappa_VDW=kappa_diff*(4.1868/10^(-2))+kappa_EU;
end    


end

