
global visc

if strcmp(visc,'Sutherland')
    PR_PROP.dat = load('CO2_80_Sutherland.dat');
elseif strcmp(visc,'JST')
    PR_PROP.dat = load(strcat('prop/','PR_CO2_',num2str(BL_PARA.p_inf,4),'_JST.dat')); 
elseif strcmp(visc,'Chung')
    PR_PROP.dat = load(strcat('prop/','PR_CO2_',num2str(BL_PARA.p_inf,4),'_redu_Chung.dat'));
elseif strcmp(visc,'constant')
    PR_PROP.dat = load('CO2_80_VDW_constant.dat');
end

PR_PROP.Zc=0.3112;
PR_PROP.R=1/PR_PROP.Zc;
PR_PROP.Ru=8.31451;
PR_PROP.M=0.0440098;
PR_PROP.omega_ac=0.224;
PR_PROP.a=0.45724;
PR_PROP.b=0.07780;
PR_PROP.dof=9;
PR_PROP.R_inf=PR_PROP.Ru./PR_PROP.M;             % Specific gas constant
PR_PROP.T_crit=3.041282000000000e+02;
PR_PROP.p_crit=7377300;
PR_PROP.Zc_RP=0.274586376;
PR_PROP.v_crit=PR_PROP.Zc_RP*PR_PROP.R_inf*PR_PROP.T_crit/PR_PROP.p_crit;

PR_PROP.h0_approx=PR_PROP.dat(:,1);
PR_PROP.T0_approx=PR_PROP.dat(:,2);
PR_PROP.p0_approx=PR_PROP.dat(:,3);
PR_PROP.rho0_approx=PR_PROP.dat(:,4);
PR_PROP.Mu_approx=PR_PROP.dat(:,5);
PR_PROP.Cp_approx=PR_PROP.dat(:,6);
PR_PROP.Kappa_approx=PR_PROP.dat(:,7);
PR_PROP.a_sound_approx=PR_PROP.dat(:,8);
PR_PROP.e0_approx=PR_PROP.dat(:,9);  

