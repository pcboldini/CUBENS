function [h_PR,rho,Mu,Cp_rho,Cv_rho,Kappa,a_sound,e_PR,PR_PROP] = PR(T,PR_PROP,BL_PARA)

global visc

CvR=BL_PARA.CvR;
a=PR_PROP.a;
b=PR_PROP.b;
omega_ac=PR_PROP.omega_ac;
Zc=PR_PROP.Zc;
Zc_1=PR_PROP.Zc^(-1);
Zc_2=PR_PROP.Zc^(-2);
p_inf=BL_PARA.p_inf;

K=0.37464+1.54226*omega_ac-0.26992*omega_ac^2;
PR_PROP.K=K;

% Density/Volume

PR_PROP.rho0=interp1(PR_PROP.T0_approx,PR_PROP.rho0_approx,T,'spline');
PR_PROP.v0=1/PR_PROP.rho0;

options = optimoptions('fsolve','Display','none');
v=fsolve(@PR_calc_v,PR_PROP.v0,options);

rho=1/v;
alpha=(1+K*(1-sqrt(T)))^2;

% Thermodynamic properties

% Specific heats

    Cv_rho=CvR-(a*K*(K+1))/(4*b*sqrt(2*T))*log((1+b*Zc_1*rho*(1-sqrt(2)))/(1+b*Zc_1*rho*(1+sqrt(2))));
    dPdT_Rho=K*sqrt(alpha/T)*(rho^2*a*Zc_2)/(1+2*b*Zc_1*rho-b^2*Zc_2*rho^2)-(rho*Zc_1)/(rho*b*Zc_1-1);
    dPdRho_T=Zc_1*T/(Zc_1*b*rho-1)^2-a*Zc_2*alpha*(2*rho+2*b*Zc_1*rho^2)/(1+2*b*Zc_1*rho-rho^2*b^2*Zc_2)^2;
    Cp_rho=Cv_rho+T/rho^2*Zc*dPdT_Rho^2/dPdRho_T;

    % Internal Energy
    
    e_PR=CvR*T*Zc_1+a*Zc_1/(2*sqrt(2)*b)*((K+1)^2-K*(K+1)*sqrt(T))*log((1+b*Zc_1*rho*(1-sqrt(2)))/(1+b*Zc_1*rho*(1+sqrt(2))));
    
    % Enthalpy
  
    h_PR=e_PR+p_inf/rho;

    % Speed of sound

    a_sound_squared=dPdRho_T+Zc*T/rho^2/Cv_rho*dPdT_Rho^2;
    a_sound=sqrt(a_sound_squared);

    
    % Dynamic viscosity & thermal conductivity
    
    if strcmp(visc,'JST')

        [Mu,Kappa] = JST(rho,T,BL_PARA,PR_PROP); %% non-dimensional

    elseif strcmp(visc,'Chung')

        [Mu,Kappa] = Chung(rho,T,PR_PROP); %% dimensional

    elseif strcmp(visc,'Constant')
       
        Mu=1;  
        Kappa=1;
       
    end 

% Subroutine VDW

    function Vol = PR_calc_v(v)

    alpha=(1+K*(1-sqrt(T)))^2;
    Vol = T*Zc_1/(v-b*Zc_1) - alpha*a*Zc_2/(v*(v+2*b*Zc_1)-b^2*Zc_2) - p_inf;
    
    end


end

