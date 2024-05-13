%% Blasius solver

function [BL_PARA,BL_MESH,BaseFlow,Pressure,Thermo,Energy]=Exec_BL(BL_MESH,BL_PARA)

global Y_Point LST

Rg=BL_PARA.Rg;
Pr_inf=BL_PARA.Pr_inf;
gamma=BL_PARA.gamma;
R_inf=BL_PARA.R_inf;
rho_inf=BL_PARA.rho_inf;
T_inf=BL_PARA.T_inf;
U_inf=BL_PARA.U_inf;
Cp=BL_PARA.cp;
Cv=BL_PARA.cv;
cp_inf=BL_PARA.cp_inf;

% Initial Profile:

Y_Point_Ex=Y_Point/2;
    
[vars_int]=IG_RK4(BL_PARA,Y_Point);

disp('BL solved!');

% Extrapolation

eta1D=vars_int(:,6);
F0=vars_int(:,1);
F1=vars_int(:,2);
F2=vars_int(:,3);
G0=vars_int(:,4);
G1=vars_int(:,5);

% Calculation of the base flow variables

IntU=F0;
U=F1;
DU=F2;
T=G0;
Rho=1./T;

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
V_ns_E=[V_ns;V_ns(Y_Point)*ones(Y_Point_Ex,1)];

% Final interpolation of non-conservative variables

Rho0=interp1(Y_E,Rho_E,y_BL,'spline');
U0=interp1(Y_E,U_E,y_BL,'spline');
T0=interp1(Y_E,T_E,y_BL,'spline');
V0_ns=interp1(Y_E,V_ns_E,y_BL,'spline');

BL_PARA.V0_ns_inf=V0_ns(end);

if strcmp(LST,'true')

    Rho0_y=Dy_BL*Rho0;
    Rho0_yy=DDy_BL*Rho0;
    U0_y=Dy_BL*U0;
    U0_yy=DDy_BL*U0;
    T0_y=Dy_BL*T0;
    T0_yy=DDy_BL*T0;

    p = Rho0*Rg.*T0;
    p_T = Rho0*Rg;
    p_TT = zeros(length(p),1);
    p_Rho = T0*Rg;
    p_RhoRho = zeros(length(p),1);
    p_RhoT = Rg*ones(length(p),1);

    [mu,mu_T,mu_TT]   = Sutherland(T0);
    mu_Rho = zeros(length(T0),1);
    mu_RhoRho = zeros(length(T0),1);
    mu_RhoT = zeros(length(T0),1);

    kappa = mu;
    kappa_T = mu_T;
    kappa_TT = mu_TT;
    kappa_Rho = zeros(length(T0),1);
    kappa_RhoRho = zeros(length(T0),1);
    kappa_RhoT = zeros(length(T0),1);

    e = BL_PARA.cv*T0;
    e_T = BL_PARA.cv*ones(length(e),1);
    e_Rho = zeros(length(e),1);  

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
BaseFlow(:,10)=Pr_inf;
BaseFlow(:,11)=sqrt(gamma*R_inf*T0*T_inf);
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
Thermo(:,7)=kappa;
Thermo(:,8)=kappa_T;
Thermo(:,9)=kappa_TT;
Thermo(:,10)=kappa_Rho;
Thermo(:,11)=kappa_RhoRho;
Thermo(:,12)=kappa_RhoT;

Energy=zeros(length(y_BL),3); % [e,e_T,e_rho]
Energy(:,1)=e;
Energy(:,2)=e_T;
Energy(:,3)=e_Rho;

else

BaseFlow=zeros(length(y_BL),6); % [rho,u,T]
BaseFlow(:,1)=Rho0;
BaseFlow(:,2)=U0;
BaseFlow(:,3)=T0;
BaseFlow(:,4)=Pr_inf;
BaseFlow(:,5)=sqrt(gamma*R_inf*T0*T_inf)/U_inf;
BaseFlow(:,6)=V0_ns;

Pressure=zeros(length(y_BL),1); % [p]
Pressure(:,1)=Rho0.*Rg.*T0;

Thermo=zeros(length(y_BL),2); % [mu,kappa]
for i=1:length(y_BL)
[mu(i),~,~]   = Sutherland(T0(i));
kappa(i)=mu(i);
end

Thermo(:,1)=mu;
Thermo(:,2)=kappa;

Energy=zeros(length(y_BL),1); % [e]
Energy(:,1)=Cv*T0;

% Isobaric expansion coefficient
alpha_v=1./T0;
fact_BC=Cp./alpha_v;

end

end


