function [BL_PARA,VDW_PROP]=Set_BL_Para_IG

global wall_bc calc_Ec visc

% Free stream values
if strcmp(calc_Ec,'true')
    BL_PARA.M_inf     = 0.1;               % Mach number
elseif strcmp(calc_Ec,'fixed')
    BL_PARA.Ec_inf     = 0.05;                % Eckert number
end
BL_PARA.p_inf     = 1000;   % Pressure (Pa)
BL_PARA.T_inf     = 300; % Free stream temperature (K) 
BL_PARA.Pr_inf    = 0.75; % Prandtl number (Pr)
BL_PARA.dof       = 9;
BL_PARA.CvR       = BL_PARA.dof/2;    

visc='PowerLaw'; % Constant/PowerLaw/Sutherland

% Wall boundary condition
wall_bc='adiab'; % isoth_nrbc
if strcmp(wall_bc,'adiab')
    BL_PARA.T_wall    = -1.0;                  % Wall temperature
    else if strcmp(wall_bc,'isoth')
    BL_PARA.T_wall    = 345;  %1.0523       % Wall temperature
    end
end

%% IG

BL_PARA.gamma     = 1.4;                % Specific heat ratio
BL_PARA.R_inf     = 287.8;                % Specific gas constant

BL_PARA.cp_inf    = BL_PARA.R_inf*BL_PARA.gamma/(BL_PARA.gamma-1);
BL_PARA.a_inf     = sqrt(BL_PARA.R_inf*BL_PARA.gamma*BL_PARA.T_inf);

% Eckert number
if strcmp(calc_Ec,'true')
    BL_PARA.Ec_inf=(BL_PARA.gamma -1)*BL_PARA.M_inf^2;              % Eckert number
elseif strcmp(calc_Ec,'fixed')
    BL_PARA.M_inf     = sqrt( BL_PARA.Ec_inf/(BL_PARA.gamma -1) );                  % Mach number    
end

% Free stream calculation

BL_PARA.rho_inf   = BL_PARA.p_inf/BL_PARA.R_inf/BL_PARA.T_inf; % Density
BL_PARA.h_inf   = BL_PARA.cp_inf*BL_PARA.T_inf;
BL_PARA.Rg        = 1/BL_PARA.gamma/BL_PARA.M_inf^2;
BL_PARA.cp        = BL_PARA.Rg*BL_PARA.gamma/(BL_PARA.gamma-1);                 % specific heat capacity (const. p)
BL_PARA.cv        = BL_PARA.cp/BL_PARA.gamma;             % specific heat capacity (const. V)

BL_PARA.U_inf=BL_PARA.M_inf * sqrt(BL_PARA.gamma * BL_PARA.R_inf * BL_PARA.T_inf);   % Free stream velocity
BL_PARA.Z_inf  = BL_PARA.p_inf/BL_PARA.rho_inf/BL_PARA.R_inf/BL_PARA.T_inf;
% Viscosity (Sutherland's law)

BL_PARA.T_ref     = 273; % Sutherland's law constant
BL_PARA.S_mu_ref         = 111; % Sutherland's law constant
BL_PARA.mu_ref    = 1.716*10^(-5); % Sutherland's law reference viscosity (Pa s)

BL_PARA.kappa_ref    = 0.0241; % Sutherland's law reference viscosity (Pa s)
BL_PARA.S_kappa_ref        = 194; % Sutherland's law constant

BL_PARA.powerexp=0.75;

[BL_PARA.mu_inf,BL_PARA.kappa_inf]   = Sutherland_dim(BL_PARA);

BL_PARA.Re_unit   = BL_PARA.U_inf*BL_PARA.rho_inf/BL_PARA.mu_inf; % Unit Reynolds number

VDW_PROP=1;

end

