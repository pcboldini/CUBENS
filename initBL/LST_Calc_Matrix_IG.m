    function [Lt,Lx,Ly,Lz,Lq,Vxx,Vxy,Vyy,Vxz,Vyz,Vzz]=LST_Calc_Matrix_IG(Re0,Pr,Ec,M,FLOW_IG_LST,PRESSURE_IG_LST,THERMO_IG_LST,ENERGY_IG_LST)

    % Parallel Assumption
    
    % Flow variables
    
    Rho=FLOW_IG_LST(1);
    Rho_y=FLOW_IG_LST(2);
    Rho_yy=FLOW_IG_LST(3);
    U=FLOW_IG_LST(4);
    U_y=FLOW_IG_LST(5);
    U_yy=FLOW_IG_LST(6);
    T=FLOW_IG_LST(7);
    T_y=FLOW_IG_LST(8);
    T_yy=FLOW_IG_LST(9);
    
    p=PRESSURE_IG_LST(1);
    p_T=PRESSURE_IG_LST(2);
    p_TT=PRESSURE_IG_LST(3);
    p_rho=PRESSURE_IG_LST(4);
    p_rhorho=PRESSURE_IG_LST(5);
    p_rhoT=PRESSURE_IG_LST(6);
    
    % Thermodynamic properties
    
    mu=THERMO_IG_LST(1);
    mu_T=THERMO_IG_LST(2);
	mu_TT=THERMO_IG_LST(3);
	mu_rho=THERMO_IG_LST(4);
	mu_rhorho=THERMO_IG_LST(5);
	mu_rhoT=THERMO_IG_LST(6);
	kappa=THERMO_IG_LST(7);
	kappa_T=THERMO_IG_LST(8);
	kappa_TT=THERMO_IG_LST(9);
	kappa_rho=THERMO_IG_LST(10);
	kappa_rhorho=THERMO_IG_LST(11);
	kappa_rhoT=THERMO_IG_LST(12);
    
    mu_y = mu_T * T_y + mu_rho* Rho_y;
    kappa_y = kappa_T * T_y + kappa_rho* Rho_y;
    
   % Internal energy
   
    e=ENERGY_IG_LST(1);
    e_T=ENERGY_IG_LST(2);
	e_rho=ENERGY_IG_LST(3);
    
    e_y = e_T*T_y + e_rho*Rho_y; 

    factp=1;
    facte=1;
    factK=Ec;

    %-----------------------------------Lt--------------------------------------
    Lt = diag([1, Rho, Rho, Rho, Rho*e_T]);
    Lt(5,1) = Rho*e_rho;
    %-----------------------------------Lx--------------------------------------
    Lx = zeros(5,5);
    % Continuum
    Lx(1,1) = U;
    Lx(1,2) = Rho;
    % x Momentum
    Lx(2,1) = p_rho/factp;
    Lx(2,2) = Rho*U;
    Lx(2,3) =-mu_y/Re0;
    Lx(2,5) = p_T/factp;
    %  y Momentum
    Lx(3,1) = -mu_rho*U_y/Re0;
    Lx(3,2) = 2/3*mu_y/Re0;
    Lx(3,3) = Rho*U;
    Lx(3,5) =-mu_T*U_y/Re0;
    %  z Momentum
    Lx(4,4) = Rho*U;
    % Energy
    Lx(5,1) = Rho*U*e_rho;
    Lx(5,2) = p*facte/factp;
    Lx(5,3) =-2*mu*facte/Re0*U_y;
    Lx(5,5) = Rho*U*e_T;
    %-----------------------------------Ly--------------------------------------
    Ly = zeros(5,5);
    % Continuum
    Ly(1,3) = Rho;
    % x Momentum
    Ly(2,1) = -mu_rho/Re0*U_y;
    Ly(2,2) = -mu_y/Re0;
    Ly(2,5) = -mu_T/Re0*U_y;
    %  y Momentum
    Ly(3,1) = p_rho/factp;
    Ly(3,3) = -2*mu_y/Re0+2/3*mu_y/Re0;
    Ly(3,5) = p_T/factp;
    %  z Momentum
    Ly(4,4) = -mu_y/Re0;
    % Energy
    Ly(5,1) = -kappa_rho*T_y/Re0/Pr/factK;
    Ly(5,2) =-2*mu*facte/Re0*U_y;
    Ly(5,3) = facte*p/factp;
    Ly(5,5) = -(kappa_y+kappa_T*T_y)/Re0/Pr/factK;
    %-----------------------------------Lz--------------------------------------
    Lz = zeros(5,5);
    %Continuum
    Lz(1,4) = Rho;
    % x Momentum
    % y Momentum
    Lz(3,4) = 2/3*mu_y/Re0;
    % z Momentum
    Lz(4,1) = p_rho/factp;
    Lz(4,3) =-mu_y/Re0;
    Lz(4,5) = p_T/factp;
    % Energy
    Lz(5,4) = facte*p/factp;
    %-----------------------------------Lq--------------------------------------
    Lq = zeros(5,5);
    % Continuum
    Lq(1,3) = Rho_y;
    %  x Momentum
    Lq(2,1) = -mu_rho*U_yy/Re0 - U_y/Re0*( mu_rhorho*Rho_y+mu_rhoT*T_y ) ;% 
    Lq(2,3) = Rho*U_y;
    Lq(2,5) = - 1/Re0*mu_T*U_yy - 1/Re0*U_y*( mu_TT*T_y + mu_rhoT*Rho_y );
    %  y Momentum
    Lq(3,1) = p_rhorho/factp*Rho_y+p_rhoT/factp*T_y;
    Lq(3,5) = p_TT/factp*T_y+p_rhoT/factp*Rho_y;
    % z Momentum
    % Energy
    Lq(5,1) = -1/Re0/Pr/factK*( T_yy*kappa_rho + T_y*( kappa_rhorho*Rho_y + kappa_rhoT*T_y ) )...
             -1/Re0*mu_rho*U_y^2*facte;
    Lq(5,3) = Rho*e_y;
    Lq(5,5) = -1/Re0/Pr/factK*( T_yy*kappa_T + T_y*( kappa_TT*T_y + kappa_rhoT*Rho_y ) )...
             -1/Re0*mu_T*U_y^2*facte;

    %----------------------------------Vxx-------------------------------------
    Vxx = diag([0,-4/3*mu/Re0,-mu/Re0,-mu/Re0,-kappa/Re0/Pr/factK]); %
    %----------------------------------Vyy-------------------------------------
    Vyy = diag([0,-mu/Re0,-4/3*mu/Re0,-mu/Re0,-kappa/Re0/Pr/factK]);
    %----------------------------------Vzz-------------------------------------
    Vzz = diag([0,-mu/Re0,-mu/Re0,-4/3*mu/Re0,-kappa/Re0/Pr/factK]);
    %----------------------------------Vxy-------------------------------------
    Vxy = zeros(5,5);
    Vxy(2,3) = -1/3*mu/Re0;
    Vxy(3,2) = -1/3*mu/Re0;
    %----------------------------------Vxz-------------------------------------
    Vxz = zeros(5,5);
    Vxz(2,4) = -1/3*mu/Re0;
    Vxz(4,2) = -1/3*mu/Re0;
    %----------------------------------Vyz-------------------------------------
    Vyz = zeros(5,5);
    Vyz(3,4) = -1/3*mu/Re0;
    Vyz(4,3) = -1/3*mu/Re0;
    %--------------------------------------------------------------------------

    end