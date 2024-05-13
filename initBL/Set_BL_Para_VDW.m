function [BL_PARA,VDW_PROP]=Set_BL_Para_VDW

global wall_bc calc_Ec visc

% Free stream values
if strcmp(calc_Ec,'true')
    BL_PARA.M_inf     = 0.1;   %0.2               % Mach number
elseif strcmp(calc_Ec,'fixed')
    BL_PARA.Ec_inf     = 0.05;   %0.3               % Eckert number
end
BL_PARA.p_inf     = 1.1; % Reduced Pressure 0.1356
BL_PARA.T_inf     = 0.9; % Reduced Free stream temperature  0.9207
BL_PARA.Pr_inf    = 1.0; % Prandtl number  3.154466
BL_PARA.CvR       = 9/2;   %% IG  

visc='Chung'; % Constant/JST/Chung

% Wall boundary condition
wall_bc='adiab'; % isoth_nrbc
if strcmp(wall_bc,'adiab')
    BL_PARA.T_wall    = -1;                  % Wall temperature
    else if strcmp(wall_bc,'isoth')
    BL_PARA.T_wall    = 1.1;  %1.0523       % Wall temperature
    end
end

%% VDW

% Loading EoS properties
Load_VDW_Prop;

if strcmp(visc,'Chung')
    VDW_PROP=Load_Chung(VDW_PROP);
end

% Free stream calculation
[BL_PARA.h_inf,BL_PARA.Rho_inf,BL_PARA.mu_inf,BL_PARA.Cp_inf,BL_PARA.kappa_inf,BL_PARA.a_inf,BL_PARA.e_inf,VDW_PROP]   = VDW(BL_PARA.T_inf,VDW_PROP,BL_PARA);

% Free stream calculation
BL_PARA.gamma = 1.3;
BL_PARA.Z_inf  = BL_PARA.p_inf/BL_PARA.Rho_inf/VDW_PROP.R/BL_PARA.T_inf;

if strcmp(calc_Ec,'true')
BL_PARA.Ec_inf=BL_PARA.M_inf^2*BL_PARA.a_inf^2/BL_PARA.Cp_inf/BL_PARA.T_inf*VDW_PROP.Zc; % Mach number
elseif strcmp(calc_Ec,'fixed')
BL_PARA.M_inf=sqrt(BL_PARA.Ec_inf*BL_PARA.Cp_inf*BL_PARA.T_inf/BL_PARA.a_inf^2/VDW_PROP.Zc); % Mach number
end

[BL_PARA.Cp_max,Cp_max_index]=max(VDW_PROP.Cp_approx);
BL_PARA.T_pc=VDW_PROP.T0_approx(Cp_max_index);

[~,BL_PARA.Rho_pc,~,~,~,~,BL_PARA.e_pc] = VDW(BL_PARA.T_pc,VDW_PROP,BL_PARA);

BL_PARA.Cp=BL_PARA.Cp_inf*BL_PARA.T_inf/BL_PARA.M_inf^2/BL_PARA.a_inf^2/VDW_PROP.Zc;

end

