%% COMPRESSIBLE BLASIUS SOLVER

    %  - Validated with Fedorov (FEDEROV,A.& TUMIN,A. High-speed boundary-layer
    %  instability: old terminology and a newframework.AIAA J.49(8), 1647–1657,
    %  2011): Ma=4.2, T_0=300K, Re=2000, adiabatic
    
    %  - Validated with Ma (MA, Y. & ZHONG, X. Receptivity of a supersonic 
    % boundary layer over a flat plate. Part 1. Wave structures and interactions.
    % J. Fluid Mech. 488, 31–78, 2003): Ma=4.5, T_0=329K, Re=2000, adiabatic
    
clc
clear all
close all

%% Parameters 

global save_output show_plot save_plot 
global USE_EOS fluid wall_bc calc_Ec visc
global Y_Point LST

USE_EOS='PR'; % IG/VdW/PR

save_output='true'; % true/false
show_plot='false'; % true/false
save_plot='false'; % true/false
LST="false";

fluid='-1';
% Keep fixed Ma or Ec
calc_Ec='true'; % true/fixed

if strcmp(USE_EOS,'IG')
    [BL_PARA,VDW_PROP]=Set_BL_Para_IG;
elseif strcmp(USE_EOS,'VdW')
    [BL_PARA,VDW_PROP]=Set_BL_Para_VDW;
elseif strcmp(USE_EOS,'PR')
    [BL_PARA,PR_PROP]=Set_BL_Para_PR;
end

%% Mesh

% Grid geometry input
BL_MESH.N_BL= 2000;         % number of collocation points in y
BL_MESH.ymax_factor = 50;   % domain size in y-times the boundary layer thickness

%% Call solver 
Y_Point = 20000;

if strcmp(USE_EOS,'IG')
[BL_PARA,BL_MESH,FLOW,PRESSURE,THERMO,ENERGY]=Exec_BL_IG(BL_MESH,BL_PARA);
elseif strcmp(USE_EOS,'VdW')
[FLOW,PRESSURE,THERMO,ENERGY,BL_PARA,BL_MESH]=Exec_BL_VDW(BL_MESH,BL_PARA,VDW_PROP);
elseif strcmp(USE_EOS,'PR')
[FLOW,PRESSURE,THERMO,ENERGY,BL_PARA,BL_MESH]=Exec_BL_PR(BL_MESH,BL_PARA,PR_PROP);
end

if strcmp(LST,'true')
    LST_PARA=Set_LST_Parameters;

    LST_MESH.N_LST= 301;         % number of collocation points in y
    LST_MESH.ymax_LST  = 300;                   % domain size in y as a multiple of the incompressible boundary layer thickness
    LST_MESH.yi_LST     = 10;                   % grid clustering parameter

    LST_MESH=Set_LST_Grid_new(LST_MESH);
    Ny_ext =2000;
    if strcmp(USE_EOS,'IG')
        [FLOW_LST,PRESSURE_LST,THERMO_LST,ENERGY_LST] = BaseFlowVariables_Interpol_IG(BL_MESH.N_BL,Ny_ext,BL_MESH.y_BL,LST_MESH.y_LST,BL_PARA,FLOW,PRESSURE,THERMO,ENERGY);
        [LST_SPECTRUM] = LST_Calc_IG(BL_PARA,LST_PARA,LST_MESH,FLOW_LST,PRESSURE_LST,THERMO_LST,ENERGY_LST);
    elseif strcmp(USE_EOS,'VdW')
        [FLOW_LST,PRESSURE_LST,THERMO_LST,ENERGY_LST] = BaseFlowVariables_Interpol_VdW(BL_MESH.N_BL,Ny_ext,BL_MESH.y_BL,LST_MESH.y_LST,BL_PARA,FLOW,PRESSURE,THERMO,ENERGY);
        [LST_SPECTRUM] = LST_Calc_VdW(BL_PARA,LST_PARA,LST_MESH,FLOW_LST,PRESSURE_LST,THERMO_LST,ENERGY_LST);
    elseif strcmp(USE_EOS,'PR')
        [FLOW_LST,PRESSURE_LST,THERMO_LST,ENERGY_LST] = BaseFlowVariables_Interpol_PR(BL_MESH.N_BL,Ny_ext,BL_MESH.y_BL,LST_MESH.y_LST,BL_PARA,FLOW,PRESSURE,THERMO,ENERGY);
        [LST_SPECTRUM] = LST_Calc_PR(BL_PARA,LST_PARA,LST_MESH,FLOW_LST,PRESSURE_LST,THERMO_LST,ENERGY_LST);
   
    end
    [LST_EIGS]=Postproc_LST_Spectrum(BL_PARA,LST_PARA,LST_MESH,LST_SPECTRUM,FLOW_LST,PRESSURE_LST,THERMO_LST,ENERGY_LST);
end

%% Postprocessing DNS
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('DNS parameters:');
fprintf ('USE_EOS=%.3s \n',USE_EOS);
fprintf ('Ec= %.3e \n',BL_PARA.Ec_inf);
fprintf ('Pra= %.3f \n',BL_PARA.Pr_inf);
fprintf ('gam= %.3f \n',BL_PARA.gamma);
fprintf ('Ma= %.3e \n',BL_PARA.M_inf);
if strcmp(USE_EOS,'IG')
    fprintf ('Rgas= %.3f \n',BL_PARA.Rg);
    fprintf ('Pref= %.3f \n',BL_PARA.Rg);
elseif strcmp(USE_EOS,'VdW') || strcmp(USE_EOS,'PR')
    fprintf ('Pref= %.3f \n',BL_PARA.p_inf);
    fprintf ('Cv/R= %.3f \n',BL_PARA.CvR);
end
fprintf ('USE_VISC=%.10s \n\n',visc);
fprintf ('BC_wall = %s\n',wall_bc);
if strcmp(USE_EOS,'IG')
    fprintf ('Twall (K)=%.5f \n\n',BL_PARA.T_wall);
elseif strcmp(USE_EOS,'VdW') || strcmp(USE_EOS,'PR')
    fprintf ('Twall (T/Tc)=%.5f \n\n',BL_PARA.T_wall);
end

BL_MESH.y_DNS=BL_MESH.y_BL./BL_PARA.delta_eta_BL;


if strcmp(LST,'true')
    
LST_MESH.y_DNS=LST_MESH.y_LST./BL_PARA.delta_eta_BL;

% BF
x_Postproc=BL_MESH.y_DNS;
Rho_Postproc=FLOW(:,1);
T_Postproc=FLOW(:,3);
W_Postproc=FLOW(:,4);
U_Postproc=FLOW(:,12);

% LST
x_LST_Postproc=LST_MESH.y_DNS;
W_LST_Postproc=LST_EIGS.u_EF;
U_LST_Postproc=LST_EIGS.v_EF;
T_LST_Postproc=LST_EIGS.T_EF;

else

% Wall-normal direction
x_Postproc=BL_MESH.y_DNS;
% Density
Rho_Postproc=FLOW(:,1);
% Streamwise velocity
W_Postproc=FLOW(:,2);
% Wall-normal velocity
U_Postproc=FLOW(:,6);

end

fprintf ('delta_99 (start)= %.3f \n',BL_PARA.delta_eta_BL);

% Write profiles 

if strcmp(save_output,'true')

    if strcmp(LST,'true')
        saveFile_x_LST         = fopen(strcat('./inputDNS/prof_x_LST.bin'),'w');
        saveFile_T_LST         = fopen(strcat('./inputDNS/prof_T_LST.bin'),'w');
        saveFile_w_LST         = fopen(strcat('./inputDNS/prof_w_LST.bin'),'w'); % Streamwise
        saveFile_u_LST         = fopen(strcat('./inputDNS/prof_u_LST.bin'),'w'); % Wall-normal
        
        fwrite(saveFile_x_LST,x_LST_Postproc,'double');
        fwrite(saveFile_T_LST,[real(T_LST_Postproc) imag(T_LST_Postproc)],'double');
        fwrite(saveFile_w_LST,[real(W_LST_Postproc) imag(W_LST_Postproc)],'double');
        fwrite(saveFile_u_LST,[real(U_LST_Postproc) imag(U_LST_Postproc)],'double');
    end
    
saveFile_x         = fopen(strcat('./inputDNS/prof_x.bin'),'w');
saveFile_r         = fopen(strcat('./inputDNS/prof_r.bin'),'w');
saveFile_w         = fopen(strcat('./inputDNS/prof_w.bin'),'w'); % Streamwise
saveFile_u         = fopen(strcat('./inputDNS/prof_u.bin'),'w'); % Wall-normal

fwrite(saveFile_x,x_Postproc,'double');
fwrite(saveFile_r,Rho_Postproc,'double');
fwrite(saveFile_w,W_Postproc,'double');
fwrite(saveFile_u,U_Postproc,'double');

% Write BASEFLOW parameters

saveFile_h         = fopen(strcat('./inputDNS/initBL_params.h'),'w');
fprintf(saveFile_h,'INIT  = "initBoundaryLayer"\n');
fprintf(saveFile_h,'! ------------ use equation of state ------------------ \n');
fprintf(saveFile_h,strcat('USE_EOS  = "',USE_EOS,'"\n'));
fprintf(saveFile_h,'! ------------ set non-dimensional free-stream values for computation ------------------ \n');
fprintf(saveFile_h,strcat('delta99 = ', num2str(BL_PARA.delta_eta_BL,'%.10f'),'\n')); 
fprintf(saveFile_h,strcat('Pra = ', num2str(BL_PARA.Pr_inf,'%.10f'),'\n'));
fprintf(saveFile_h,strcat('Ec = ', num2str(BL_PARA.Ec_inf,'%.10e'),'\n'));
fprintf(saveFile_h,strcat('Ma = ', num2str(BL_PARA.M_inf,'%.10e'),'\n'));
if strcmp(USE_EOS,'IG')
    fprintf(saveFile_h,strcat('Tinf = ', num2str(BL_PARA.T_inf,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('Pref = ', num2str(BL_PARA.Rg,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('ig_gam = ', num2str(BL_PARA.gamma,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('eos_Rgas = ', num2str(BL_PARA.Rg,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('eos_dof = ', num2str(BL_PARA.dof,'%.6f'),'\n'));
elseif strcmp(USE_EOS,'VdW')
    fprintf(saveFile_h,strcat('eos_dof = ', num2str(VDW_PROP.dof,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('eos_ac = ', num2str(VDW_PROP.omega_ac,'%.6f'),'\n'));
    fprintf(saveFile_h,'! ------------ set dimensional free-stream values for computation ------------------ \n');
    fprintf(saveFile_h,strcat('Tcrit = ', num2str(VDW_PROP.T_crit,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('Pcrit = ', num2str(VDW_PROP.p_crit,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('Vcrit = ', num2str(VDW_PROP.v_crit,'%.6e'),'\n'));
    fprintf(saveFile_h,strcat('eos_Rgas = ', num2str(VDW_PROP.R_inf,'%.6f'),'\n'));
elseif strcmp(USE_EOS,'PR')
    fprintf(saveFile_h,strcat('eos_dof = ', num2str(PR_PROP.dof,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('eos_ac = ', num2str(PR_PROP.omega_ac,'%.6f'),'\n'));
    fprintf(saveFile_h,'! ------------ set dimensional free-stream values for computation ------------------ \n');
    fprintf(saveFile_h,strcat('Tcrit = ', num2str(PR_PROP.T_crit,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('Pcrit = ', num2str(PR_PROP.p_crit,'%.6f'),'\n'));
    fprintf(saveFile_h,strcat('Vcrit = ', num2str(PR_PROP.v_crit,'%.6e'),'\n')); %% from REFPROP for Chung
    fprintf(saveFile_h,strcat('eos_Rgas = ', num2str(PR_PROP.R_inf,'%.6f'),'\n'));
end
fprintf(saveFile_h,'! ------------ set wall BC for computation ------------------ \n');
fprintf(saveFile_h,strcat('wall_bc  = "',wall_bc,'"\n'));
if strcmp(wall_bc,'adiab')
    fprintf(saveFile_h,strcat('Twall = ', num2str(BL_PARA.T_wall,'%.6f'),'\n'));
elseif strcmp(wall_bc,'isoth')
    fprintf(saveFile_h,strcat('Twall = ', num2str(BL_PARA.T_wall./BL_PARA.T_inf,'%.6f'),'\n'));
end

if strcmp(USE_EOS,'VdW') || strcmp(USE_EOS,'PR')
fprintf(saveFile_h,'! ------------ set reference free-stream values for computation ------------------ \n');
fprintf(saveFile_h,strcat('Tref = ', num2str(BL_PARA.T_inf,'%.10f'),'\n'));
fprintf(saveFile_h,strcat('Pref = ', num2str(BL_PARA.p_inf,'%.10f'),'\n')); 
fprintf(saveFile_h,strcat('Rhoref = ', num2str(BL_PARA.Rho_inf,'%.10f'),'\n'));

fprintf(saveFile_h,strcat('Cpref = ', num2str(BL_PARA.Cp_inf,'%.10f'),'\n'));
fprintf(saveFile_h,strcat('SOSref = ', num2str(BL_PARA.a_inf,'%.10f'),'\n'));
end

if strcmp(USE_EOS,'VdW')
fprintf(saveFile_h,strcat('Rref = ', num2str(VDW_PROP.R,'%.6f'),'\n'));
elseif strcmp(USE_EOS,'PR')
fprintf(saveFile_h,strcat('Rref = ', num2str(PR_PROP.R,'%.6f'),'\n'));
end

fprintf(saveFile_h,'! ------------ set viscosity and conductivity ------------------ \n');
fprintf(saveFile_h,strcat('USE_VISC  = "',visc,'"\n'));
if strcmp(USE_EOS,'IG')
    fprintf(saveFile_h,strcat('Stref = ', num2str(BL_PARA.T_ref,'%.10f'),'\n'));
    fprintf(saveFile_h,strcat('Muinf = ', num2str(BL_PARA.mu_inf,'%.10e'),'\n'));
    fprintf(saveFile_h,strcat('Muref = ', num2str(BL_PARA.mu_ref,'%.10e'),'\n'));
    fprintf(saveFile_h,strcat('Smuref = ', num2str(BL_PARA.S_mu_ref,'%.10f'),'\n'));
    fprintf(saveFile_h,strcat('Kinf = ', num2str(BL_PARA.kappa_inf,'%.10e'),'\n'));
    fprintf(saveFile_h,strcat('Kref = ', num2str(BL_PARA.kappa_ref,'%.10e'),'\n'));
    fprintf(saveFile_h,strcat('Skref = ', num2str(BL_PARA.S_kappa_ref,'%.10f'),'\n'));
elseif strcmp(USE_EOS,'VdW') || strcmp(USE_EOS,'PR')
    fprintf(saveFile_h,strcat('Muref = ', num2str(BL_PARA.mu_inf,'%.10e'),'\n'));
    fprintf(saveFile_h,strcat('Kref = ', num2str(BL_PARA.kappa_inf,'%.10e'),'\n'));
end

save('../postproc/initBL_param.mat','BL_PARA','BL_MESH','FLOW','THERMO','ENERGY','PRESSURE')
disp('DNS parameters saved in /inputDNS and in /postproc');

else
end

%% Plotting
if strcmp(show_plot,'true')
PostProc;
else 
end

