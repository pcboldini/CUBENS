function [BL_PARA,PR_PROP,RG_PROP]=Set_BL_Para_PR

global wall_bc visc calc_Ec

% Free stream values
if strcmp(calc_Ec,'true')
    BL_PARA.M_inf     = 0.1;
elseif strcmp(calc_Ec,'fixed')
    BL_PARA.Ec_inf     = 0.01;             % Eckert Number
end

BL_PARA.p_inf     = 1.0844; % Reduced Pressure 0.1356
BL_PARA.T_inf     = 0.9207; % Reduced Free stream temperature  0.9207

visc='Chung'; % Constant/JST/Chung

% Wall boundary condition
wall_bc='adiab'; % isoth_nrbc
if strcmp(wall_bc,'adiab')
    BL_PARA.T_wall    = -1;                  % Wall temperature
    else if strcmp(wall_bc,'isoth')
    BL_PARA.T_wall    = 1.035747;  %1.0523       % Wall temperature
    end
end

%% VDW

% Loading EoS properties
Load_PR_Prop;

if strcmp(visc,'Chung')
    PR_PROP=Load_Chung(PR_PROP);
end

BL_PARA.CvR       = PR_PROP.dof/2;   %% IG  

% Free stream calculation
[BL_PARA.h_inf,BL_PARA.Rho_inf,BL_PARA.mu_inf,BL_PARA.Cp_inf,BL_PARA.Cv_inf,BL_PARA.kappa_inf,BL_PARA.a_inf,BL_PARA.e_inf,PR_PROP]   = PR(BL_PARA.T_inf,PR_PROP,BL_PARA);

BL_PARA.Z_inf  = BL_PARA.p_inf/BL_PARA.Rho_inf/PR_PROP.R/BL_PARA.T_inf;
BL_PARA.Pr_inf=BL_PARA.Cp_inf*BL_PARA.mu_inf/BL_PARA.kappa_inf*PR_PROP.R_inf; % Prandtl number
BL_PARA.gamma = 1.3;

if strcmp(calc_Ec,'true')
BL_PARA.Ec_inf=BL_PARA.M_inf^2*BL_PARA.a_inf^2/BL_PARA.Cp_inf/BL_PARA.T_inf*PR_PROP.Zc; % Mach number
elseif strcmp(calc_Ec,'fixed')
BL_PARA.M_inf=sqrt(BL_PARA.Ec_inf*BL_PARA.Cp_inf*BL_PARA.T_inf/BL_PARA.a_inf^2/PR_PROP.Zc); % Mach number
end

[BL_PARA.Cp_max,Cp_max_index]=max(PR_PROP.Cp_approx);
BL_PARA.T_pc=PR_PROP.T0_approx(Cp_max_index);

[~,BL_PARA.Rho_pc,~,~,~,~,BL_PARA.e_pc] = PR(BL_PARA.T_pc,PR_PROP,BL_PARA);

end

