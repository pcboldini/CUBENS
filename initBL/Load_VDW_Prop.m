
global visc

if strcmp(visc,'Sutherland')
    VDW_PROP.dat = load('CO2_80_Sutherland.dat');
elseif strcmp(visc,'JST') || strcmp(visc,'Chung')
    VDW_PROP.dat = load(strcat('prop/','VDW_CO2_',num2str(BL_PARA.p_inf,3),'_redu_JST.dat')); % 80
elseif strcmp(visc,'constant')
    VDW_PROP.dat = load('VDW_CO2_80_redu.dat');
end

VDW_PROP.Zc=3/8;
VDW_PROP.R=1/VDW_PROP.Zc;
VDW_PROP.a=3;
VDW_PROP.b=1/3;
VDW_PROP.Ru=8.31451;
VDW_PROP.omega_ac=0.224;
VDW_PROP.dof=9;
VDW_PROP.M=0.0440098;
VDW_PROP.R_inf=VDW_PROP.Ru./VDW_PROP.M;             % Specific gas constant
VDW_PROP.T_crit=3.041282000000000e+02;
VDW_PROP.p_crit=7377300;
VDW_PROP.v_crit=VDW_PROP.Zc*VDW_PROP.R_inf*VDW_PROP.T_crit/VDW_PROP.p_crit;

VDW_PROP.h0_approx=VDW_PROP.dat(:,1);
VDW_PROP.T0_approx=VDW_PROP.dat(:,2);
VDW_PROP.p0_approx=VDW_PROP.dat(:,3);
VDW_PROP.rho0_approx=VDW_PROP.dat(:,4);
VDW_PROP.Mu_approx=VDW_PROP.dat(:,5);
VDW_PROP.Cp_approx=VDW_PROP.dat(:,6);
VDW_PROP.Kappa_approx=VDW_PROP.dat(:,7);
VDW_PROP.a_sound_approx=VDW_PROP.dat(:,8);
VDW_PROP.e0_approx=VDW_PROP.dat(:,9);
