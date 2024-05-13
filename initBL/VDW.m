function [h_VDW,rho,mu_VDW,Cp_VDW,kappa_VDW,a_VDW,e_VDW,VDW_PROP] = VDW(T,VDW_PROP,BL_PARA)

global visc

R=VDW_PROP.R;
a=VDW_PROP.a;
b=VDW_PROP.b;
CvR=BL_PARA.CvR;
p_inf=BL_PARA.p_inf;

VDW_PROP.rho0=interp1(VDW_PROP.T0_approx,VDW_PROP.rho0_approx,T,'spline');

options = optimoptions('fsolve','Display','none');
rho=fsolve(@VDW_calc_rho,VDW_PROP.rho0,options);

    function Rho = VDW_calc_rho(rho)
    
    Rho = (8*T)/(3*(1/rho)-1) - 3*rho^2 - p_inf ;
    
    end

sol_v=1/rho;

% Thermodynamic properties

    % Specific heats

    Cp_VDW=CvR+1/(1-(3*sol_v-1).^2/(4*T*sol_v.^3));

    % Internal energy
    
    e_VDW=R*CvR*T-a*rho;
    
    % Dynamic viscosity
    
    if strcmp(visc,'JST') 

        [mu_VDW,kappa_VDW] = JST(rho,T,BL_PARA,VDW_PROP); %% non-dimensional

    elseif strcmp(visc,'Chung') 

        [mu_VDW,kappa_VDW] = Chung(rho,T,VDW_PROP); %% dimensional
        
    elseif strcmp(visc,'Constant')
       
        mu_VDW=1;  
        kappa_VDW=1;
       
    end 
    
    % Enthalpy
    
    h_VDW=e_VDW+p_inf/rho;
    
    % Speed of sound
    
    a_sound_squared=(1+1/CvR)*R*T*(3*sol_v/(3*sol_v-1))^2-6/sol_v;
    a_VDW=sqrt(a_sound_squared);

end

