    function [Lt,Lx,Ly,Lz,Lq,Vxx,Vxy,Vyy,Vxz,Vyz,Vzz]=LST_Calc_Matrix_VdW(Re0,Pr,Ec,FLOW_VDW_LST,PRESSURE_VDW_LST,THERMO_VDW_LST,ENERGY_VDW_LST)

    % Parallel Assumption
    
    % Flow variables
    
    Rho=FLOW_VDW_LST(1);
    Rho_y=FLOW_VDW_LST(2);
    Rho_yy=FLOW_VDW_LST(3);
    U=FLOW_VDW_LST(4);
    U_y=FLOW_VDW_LST(5);
    U_yy=FLOW_VDW_LST(6);
    T=FLOW_VDW_LST(7);
    T_y=FLOW_VDW_LST(8);
    T_yy=FLOW_VDW_LST(9);
    
    p=PRESSURE_VDW_LST(1);
    p_T=PRESSURE_VDW_LST(2);
    p_TT=PRESSURE_VDW_LST(3);
    p_rho=PRESSURE_VDW_LST(4);
    p_rhorho=PRESSURE_VDW_LST(5);
    p_rhoT=PRESSURE_VDW_LST(6);
    
    % Thermodynamic properties
    
    mu=THERMO_VDW_LST(1);
    mu_T=THERMO_VDW_LST(2);
	mu_TT=THERMO_VDW_LST(3);
	mu_rho=THERMO_VDW_LST(4);
	mu_rhorho=THERMO_VDW_LST(5);
	mu_rhoT=THERMO_VDW_LST(6);
	kappa=THERMO_VDW_LST(7);
	kappa_T=THERMO_VDW_LST(8);
	kappa_TT=THERMO_VDW_LST(9);
	kappa_rho=THERMO_VDW_LST(10);
	kappa_rhorho=THERMO_VDW_LST(11);
	kappa_rhoT=THERMO_VDW_LST(12);
    
    mu_y = mu_T * T_y + mu_rho* Rho_y;
    kappa_y = kappa_T * T_y + kappa_rho* Rho_y;
    
   % Internal energy
   
    e=ENERGY_VDW_LST(1);
    e_T=ENERGY_VDW_LST(2);
	e_rho=ENERGY_VDW_LST(3);
    
    e_y = e_T*T_y + e_rho*Rho_y; 
    
    %-----------------------------------Lt--------------------------------------
    Lt = diag([1, Rho, Rho, Rho, Rho*e_T,0]);
    Lt(5,1) = Rho*e_rho;
    %-----------------------------------Lx--------------------------------------
    Lx = zeros(6,6);
    % Continuum
    Lx(1,1) = U;
    Lx(1,2) = Rho;
    % x Momentum
%    Lx(2,1) = p_rho;
    Lx(2,2) = Rho*U;
    Lx(2,3) =-mu_y/Re0;
%    Lx(2,5) = p_T;
    Lx(2,6) = 1;
    %  y Momentum
    Lx(3,1) = -mu_rho*U_y/Re0;
    Lx(3,2) = 2/3*mu_y/Re0;
    Lx(3,3) = Rho*U;
    Lx(3,5) =-mu_T*U_y/Re0;
    %  z Momentum
    Lx(4,4) = Rho*U;
    % Energy
    Lx(5,1) = Rho*U*e_rho;
    Lx(5,2) = p;
    Lx(5,3) =-2*mu/Re0*U_y;
    Lx(5,5) = Rho*U*e_T;
    %-----------------------------------Ly--------------------------------------
    Ly = zeros(6,6);
    % Continuum
    Ly(1,3) = Rho;
    % x Momentum
    Ly(2,1) = -mu_rho/Re0*U_y;
    Ly(2,2) = -mu_y/Re0;
    Ly(2,5) = -mu_T/Re0*U_y;
    %  y Momentum
%    Ly(3,1) = p_rho;
    Ly(3,3) = -2*mu_y/Re0+2/3*mu_y/Re0;
%    Ly(3,5) = p_T;
    Ly(3,6) = 1;
    %  z Momentum
    Ly(4,4) = -mu_y/Re0;
    % Energy
    Ly(5,1) = -kappa_rho*T_y/Re0/Pr/Ec;
    Ly(5,2) =-2*mu/Re0*U_y;
    Ly(5,3) = p;
    Ly(5,5) = -(kappa_y+kappa_T*T_y)/Re0/Pr/Ec;
    %-----------------------------------Lz--------------------------------------
    Lz = zeros(6,6);
    %Continuum
    Lz(1,4) = Rho;
    % x Momentum
    % y Momentum
    Lz(3,4) = 2/3*mu_y/Re0;
    % z Momentum
%    Lz(4,1) = p_rho;
    Lz(4,3) =-mu_y/Re0;
%    Lz(4,5) = p_T;
    Lz(4,6) = 1;
    % Energy
    Lz(5,4) = p;
    %-----------------------------------Lq--------------------------------------
    Lq = zeros(6,6);
    % Continuum
    Lq(1,3) = Rho_y;
    %  x Momentum
    Lq(2,1) = -mu_rho*U_yy/Re0 - U_y/Re0*( mu_rhorho*Rho_y+mu_rhoT*T_y ) ;% 
    Lq(2,3) = Rho*U_y; %%%%
    Lq(2,5) = - 1/Re0*mu_T*U_yy - 1/Re0*U_y*( mu_TT*T_y + mu_rhoT*Rho_y );
    %  y Momentum
    %Lq(3,1) = p_rhorho*Rho_y+p_rhoT*T_y;
    %Lq(3,5) = p_TT*T_y+p_rhoT*Rho_y;
    % z Momentum
    % Energy
    Lq(5,1) = -1/Re0/Pr/Ec*( T_yy*kappa_rho + T_y*( kappa_rhorho*Rho_y + kappa_rhoT*T_y ) )...
             -1/Re0*mu_rho*U_y^2;
    Lq(5,3) = Rho*e_y; %%%%
    Lq(5,5) = -1/Re0/Pr/Ec*( T_yy*kappa_T + T_y*( kappa_TT*T_y + kappa_rhoT*Rho_y ) )...
             -1/Re0*mu_T*U_y^2;
    
    % EOS
    Lq(6,1) = -p_rho;
    Lq(6,5) = -p_T;
    Lq(6,6) = 1;
         
         

    %----------------------------------Vxx-------------------------------------
    Vxx = diag([0,-4/3*mu/Re0,-mu/Re0,-mu/Re0,-kappa/Re0/Pr/Ec,0]); %
    %----------------------------------Vyy-------------------------------------
    Vyy = diag([0,-mu/Re0,-4/3*mu/Re0,-mu/Re0,-kappa/Re0/Pr/Ec,0]);
    %----------------------------------Vzz-------------------------------------
    Vzz = diag([0,-mu/Re0,-mu/Re0,-4/3*mu/Re0,-kappa/Re0/Pr/Ec,0]);
    %----------------------------------Vxy-------------------------------------
    Vxy = zeros(6,6);
    Vxy(2,3) = -1/3*mu/Re0;
    Vxy(3,2) = -1/3*mu/Re0;
    %----------------------------------Vxz-------------------------------------
    Vxz = zeros(6,6);
    Vxz(2,4) = -1/3*mu/Re0;
    Vxz(4,2) = -1/3*mu/Re0;
    %----------------------------------Vyz-------------------------------------
    Vyz = zeros(6,6);
    Vyz(3,4) = -1/3*mu/Re0;
    Vyz(4,3) = -1/3*mu/Re0;
    %--------------------------------------------------------------------------

    end