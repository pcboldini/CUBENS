%% Blasius solver

function [BaseFlow,Pressure,Thermo,Energy,BL_PARA,BL_MESH]=Exec_BL_PR(BL_MESH,BL_PARA,PR_PROP)

global save_deri wall_bc RG_model calc_deri
global visc
global Y_Point LST

Zc=PR_PROP.Zc;
Zc_1=PR_PROP.Zc^(-1);
Zc_2=PR_PROP.Zc^(-2);
a=PR_PROP.a;
b=PR_PROP.b;
K=PR_PROP.K;

T_inf=BL_PARA.T_inf;
Rho_inf=BL_PARA.Rho_inf;
Ec_inf=BL_PARA.Ec_inf;
e_inf=BL_PARA.e_inf;
Cp_inf=BL_PARA.Cp_inf;
p_inf=BL_PARA.p_inf;


%% Initial Profile
Y_Point_Ex=Y_Point/2;

[vars_int,CSol,KappaSol,PrSol,RhoSol,TSol]=PR_RK4(BL_PARA,PR_PROP,Y_Point);

disp('BL solved!');
    
% Extrapolation    

eta1D=vars_int(:,6);
F0=vars_int(:,1);
F1=vars_int(:,2);

F2=vars_int(:,3);
G0=vars_int(:,4);
G1=vars_int(:,5);

% Calculation of the base flow variables

U=F1;
T=TSol;
Rho=RhoSol;
mu=CSol./RhoSol;
kappa=KappaSol;
Pr=PrSol;
Cp=zeros(numel(T),1);
a_sound_dim=zeros(numel(T),1);
e_int=zeros(numel(T),1);

for i=1:numel(T)    
[~,~,~,Cp(i),~,a_sound_dim(i),e_int(i),~]  = PR(T(i)*T_inf,PR_PROP,BL_PARA);
end

A_sound_r=Zc/T_inf/Ec_inf/Cp_inf;
a_sound=a_sound_dim*sqrt(A_sound_r);
e_int=e_int./e_inf;

V_ns=(U.*eta1D)./(sqrt(4))- (F0)./(Rho*sqrt(2)); %
[~,index_eta_BL]=min(abs(U-0.99));

BL_PARA.delta_eta_BL=eta1D(index_eta_BL);
BL_MESH.yi_BL     = 2*BL_PARA.delta_eta_BL;              % grid clustering parameter (twice the BL)

% Calculation of the mesh based on ymax_factor times the BL thickness

BL_MESH.ymax_BL=BL_PARA.delta_eta_BL*BL_MESH.ymax_factor;

BL_MESH=Set_Grid_new(BL_MESH);
y_BL=BL_MESH.y_BL;
Dy_BL=BL_MESH.Dy_BL;
DDy_BL=BL_MESH.DDy_BL;

% Extend

y_plus=linspace(eta1D(Y_Point),BL_MESH.ymax_BL,Y_Point_Ex+1)';
Y_E=[eta1D;y_plus(2:Y_Point_Ex+1)];

U_E=[U;U(Y_Point)*ones(Y_Point_Ex,1)];
T_E=[T;T(Y_Point)*ones(Y_Point_Ex,1)];
Rho_E=[Rho;Rho(Y_Point)*ones(Y_Point_Ex,1)];
Pr_E=[Pr;Pr(Y_Point)*ones(Y_Point_Ex,1)];
Cp_E=[Cp;Cp(Y_Point)*ones(Y_Point_Ex,1)];
a_sound_E=[a_sound;a_sound(Y_Point)*ones(Y_Point_Ex,1)];
e_int_E=[e_int;e_int(Y_Point)*ones(Y_Point_Ex,1)];
V_ns_E=[V_ns;V_ns(Y_Point)*ones(Y_Point_Ex,1)];

mu_E=[mu;mu(Y_Point)*ones(Y_Point_Ex,1)];
kappa_E=[kappa;kappa(Y_Point)*ones(Y_Point_Ex,1)];

% Final interpolation of non-conservative variables

Rho0=interp1(Y_E,Rho_E,y_BL,'spline');
U0=interp1(Y_E,U_E,y_BL,'spline');
T0=interp1(Y_E,T_E,y_BL,'spline');
Pr0=interp1(Y_E,Pr_E,y_BL,'spline');
Cp0=interp1(Y_E,Cp_E,y_BL,'spline');
a_sound_0=interp1(Y_E,a_sound_E,y_BL,'spline');
e0=interp1(Y_E,e_int_E,y_BL,'spline');
V0_ns=interp1(Y_E,V_ns_E,y_BL,'spline');

mu0=interp1(Y_E,mu_E,y_BL,'spline');
kappa0=interp1(Y_E,kappa_E,y_BL,'spline');

BL_PARA.V0_ns_inf=V0_ns(end);

if strcmp(LST,'true')

    Rho0_y=Dy_BL*Rho0;
    Rho0_yy=DDy_BL*Rho0;
    U0_y=Dy_BL*U0;
    U0_yy=DDy_BL*U0;
    T0_y=Dy_BL*T0;
    T0_yy=DDy_BL*T0;

    [epsilon,p,p_T,p_TT,p_Rho,p_RhoRho,p_RhoT,mu,mu_T,mu_TT,mu_Rho,mu_RhoRho,mu_RhoT,kappa,kappa_T,kappa_TT,kappa_Rho,kappa_RhoRho,kappa_RhoT,e,e_T,e_Rho] =Calc_TP_AnaDeri_PR(BL_PARA,BL_MESH,PR_PROP,Rho0,T0,e0,mu0,kappa0);
    BL_PARA.epsilon=epsilon;

end

% Save output

if strcmp(LST,'true')

    BaseFlow=zeros(length(y_BL),12); % [rho,u,T]
    BaseFlow(:,1)=Rho0;
    BaseFlow(:,2)=Rho0_y;
    BaseFlow(:,3)=Rho0_yy;
    BaseFlow(:,4)=U0;
    BaseFlow(:,5)=U0_y;
    BaseFlow(:,6)=U0_yy;
    BaseFlow(:,7)=T0;
    BaseFlow(:,8)=T0_y;
    BaseFlow(:,9)=T0_yy;
    BaseFlow(:,10)=Pr0;
    BaseFlow(:,11)=a_sound_0;
    BaseFlow(:,12)=V0_ns;

    Pressure=zeros(length(y_BL),6); % [p,p_T,p_TT,p_rho,p_rhorho,p_rhoT]
    Pressure(:,1)=p;
    Pressure(:,2)=p_T;
    Pressure(:,3)=p_TT;
    Pressure(:,4)=p_Rho;
    Pressure(:,5)=p_RhoRho;
    Pressure(:,6)=p_RhoT;

    Thermo=zeros(length(y_BL),12); % [mu,mu_T,mu_TT,mu_rho,mu_rhorho,mu_rhoT,kappa,kappa_T,kappa_TT,kappa_rho,kappa_rhorho,kappa_rhoT]
    Thermo(:,1)=mu;
    Thermo(:,2)=mu_T;
    Thermo(:,3)=mu_TT;
    Thermo(:,4)=mu_Rho;
    Thermo(:,5)=mu_RhoRho;
    Thermo(:,6)=mu_RhoT;
    Thermo(:,7)=kappa; % kappa
    Thermo(:,8)=kappa_T; % 
    Thermo(:,9)=kappa_TT; %
    Thermo(:,10)=kappa_Rho; % 
    Thermo(:,11)=kappa_RhoRho; %
    Thermo(:,12)=kappa_RhoT; %

    Energy=zeros(length(y_BL),3); % [e,e_T,e_rho]
    Energy(:,1)=e;
    Energy(:,2)=e_T;
    Energy(:,3)=e_Rho;

    P_r=Zc/Rho_inf/T_inf/Ec_inf/Cp_inf;

else

    BaseFlow=zeros(length(y_BL),7); % [rho,u,T]
    BaseFlow(:,1)=Rho0;
    BaseFlow(:,2)=U0;
    BaseFlow(:,3)=T0;
    BaseFlow(:,4)=Pr0;
    BaseFlow(:,5)=a_sound_0;
    BaseFlow(:,6)=V0_ns;
    BaseFlow(:,7)=e0; 
    
    Pressure=zeros(length(y_BL),1); % [p]
    P_r=Zc/Rho_inf/T_inf/Ec_inf/Cp_inf;
    Pressure(:,1)=p_inf.*P_r.*ones(length(y_BL),1);

    Thermo=zeros(length(y_BL),2); % [mu,kappa]
    Thermo(:,1)=mu0;
    Thermo(:,2)=kappa0; % kappa
    
    Energy=zeros(length(y_BL),1); % [e]
    Energy(:,1)=e0;

end

% Isobaric expansion coefficient
Fac_cp=BL_PARA.T_inf/BL_PARA.M_inf^2/BL_PARA.a_inf^2/Zc;
Rho_r=Rho0.*Rho_inf;
T_r=T0.*T_inf;
alpha=(1+K.*(1-sqrt(T_r))).^2;
p_T=K.*sqrt(alpha./T_r).*(a*Rho_r.^2.*Zc_2)./(1+2*b*Zc_1*Rho_r-b^2*Zc_2*Rho_r.^2)-Rho_r*Zc_1./(b*Zc_1*Rho_r-1);
p_Rho=Zc_1*T_r./(Rho_r*b*Zc_1-1).^2-a*alpha*Zc_2.*(2*Rho_r+2*b*Zc_1*Rho_r.^2)./(1+2*b*Zc_1*Rho_r-Rho_r.^2*b^2*Zc_2).^2;

p_T=p_T.*T_inf.*P_r;
p_Rho=p_Rho.*Rho_inf.*P_r;

alpha_v=1./Rho0.*p_T./p_Rho;
fact_BC=Cp0.*Fac_cp./alpha_v;


end


