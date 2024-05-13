% Differential Matrix for Spectrum Method
 
function [epsilon,p,p_T,p_TT,p_Rho,p_RhoRho,p_RhoT,mu,mu_T,mu_TT,mu_Rho,mu_RhoRho,mu_RhoT,kappa,kappa_T,kappa_TT,kappa_Rho,kappa_RhoRho,kappa_RhoT,e,e_T,e_Rho,PT_Rho,RhoT_p] = Calc_TP_AnaDeri(BL_PARA,BL_MESH,VDW_PROP,Rho0,T0,e0,mu0,kappa0)
    global visc
    
    % Parameters
   CvR=BL_PARA.CvR;
   Zc=VDW_PROP.Zc;
   Ru=VDW_PROP.Ru;
   
   T_inf=BL_PARA.T_inf;
   Rho_inf=BL_PARA.Rho_inf;
   p_inf=BL_PARA.p_inf;
   Ec_inf=BL_PARA.Ec_inf;
   Cp_inf=BL_PARA.Cp_inf;   
   e_inf=BL_PARA.e_inf;
   mu_inf=BL_PARA.mu_inf;
   kappa_inf=BL_PARA.kappa_inf;
   kappa_factor=(0.307*CvR+0.539);
    
   N_BL=BL_MESH.N_BL;
   
   % Dimensional properties
   T_dim=T0*T_inf;
   Rho_dim=Rho0*Rho_inf;
   
   %% Initialization
   
   p_T_dim=zeros(N_BL,1);
   p_TT_dim=zeros(N_BL,1);
   p_Rho_dim=zeros(N_BL,1);
   p_RhoRho_dim=zeros(N_BL,1);
   p_RhoT_dim=zeros(N_BL,1);
   p_TRho_dim=zeros(N_BL,1);
   
   mu_dim=zeros(N_BL,1);
   mu_T_dim=zeros(N_BL,1);
   mu_TT_dim=zeros(N_BL,1);
   mu_Rho_dim=zeros(N_BL,1);
   mu_RhoRho_dim=zeros(N_BL,1);
   mu_RhoT_dim=zeros(N_BL,1);
   mu_TRho_dim=zeros(N_BL,1);
   
   kappa_dim=zeros(N_BL,1);
   kappa_T_dim=zeros(N_BL,1);
   kappa_TT_dim=zeros(N_BL,1);
   kappa_Rho_dim=zeros(N_BL,1);
   kappa_RhoRho_dim=zeros(N_BL,1);
   kappa_RhoT_dim=zeros(N_BL,1);
   kappa_TRho_dim=zeros(N_BL,1);
   
   e_dim=zeros(N_BL,1);
   e_T_dim=zeros(N_BL,1);
   e_Rho_dim=zeros(N_BL,1);

   T_Rho_dim=zeros(N_BL,1);
   T_P_dim=zeros(N_BL,1);
   
   %% Dimensional derivatives (analytical)
   
   % Pressure
   p_dim=ones(N_BL,1).*p_inf;
   p_T_dim=8*Rho_dim./(3-Rho_dim);
   p_TT_dim= zeros(length(N_BL),1);
   p_Rho_dim = 24.*T_dim./(3-Rho_dim).^2-6*Rho_dim; 
   p_RhoRho_dim = -6+48*T_dim./(3-Rho_dim).^3;
   p_RhoT_dim= 24./(3-Rho_dim).^2;
   p_TRho_dim= p_RhoT_dim;
   
   % Viscosity and thermal conductivity
  
   % Energy
   e_dim=e0.*e_inf;
   e_T_dim=CvR./Zc*ones(N_BL,1);
   e_Rho_dim=-3*ones(N_BL,1);

   % Temperature
   T_Rho_dim=1/8*(-3*p_dim./Rho_dim.^2+9-6.*Rho_dim);
   T_P_dim=1/8*(3./Rho_dim-1);
   
   %% Dimensional derivatives (numerical)
   
   T_vec=T0*T_inf;
   Rho_vec=Rho0*Rho_inf;
   
   dT = 0.01;
   dRho = 0.01;
   ddx_coeff = [-1 0 1];
   d2dx2_coeff = [1 -2 1];
   
   % Derivatives
   
   for i=1:N_BL
      T(1)= T_vec(i)-dT;
      T(2)= T_vec(i);
      T(3)= T_vec(i)+dT;
      
      Rho(1)=Rho_vec(i)-dRho;
      Rho(2)=Rho_vec(i);
      Rho(3)=Rho_vec(i)+dRho;
      
      for ii=1:length(T)
                for jj=1:length(Rho)
                if strcmp(visc,'RG')
                    
                 elseif strcmp(visc,'JST')
                     if T(ii)<=1.50
                        mu_1(ii,jj)=34*10^(-5)*T(ii)^(0.94);
                     else
                        mu_1(ii,jj)=17.78*10^(-5)*(4.58*T(ii)-1.67)^(5/8);
                     end
                     
                 f_rho(jj)=0.10230+0.023364*Rho(jj)+0.058533*Rho(jj).^2-0.040758*Rho(jj).^3+0.0093324*Rho(jj).^4;
                 mu_diff(ii,jj)=(f_rho(jj).^4-10^(-4));
                 mu_grad(ii,jj)=mu_diff(ii,jj)+mu_1(ii,jj);

                 kappa_EU(ii,jj)=15/4*Ru*mu_1(ii,jj)*kappa_factor;
   
                 kappa_1(ii,jj)=6.54*10^(-5).*(exp(0.2826.*T(ii))-1./exp(0.3976.*T(ii).^2));
                       
                 if Rho(jj)<0.50
                        f_rho(jj)=14.0*(exp(0.535*Rho(jj))-1);
                        kappa_diff(ii,jj)=(f_rho(jj)*10^(-8))/Zc^5;
                        kappa_grad(ii,jj)=kappa_diff(ii,jj)*(4.1868/10^(-2))+kappa_EU(ii,jj);
                 elseif Rho(jj)>=0.50 && Rho(jj)<2.0
                        f_rho_2(jj)=13.1*(exp(0.67*Rho(jj))-1.069);
                        kappa_diff(ii,jj)=(f_rho_2(jj)*10^(-8))/Zc^5;
                        kappa_grad(ii,jj)=kappa_diff(ii,jj)*(4.1868/10^(-2))+kappa_EU(ii,jj);
                 else
                        f_rho_3(jj)=2.976*(exp(1.155*Rho(jj))+2.016); 
                        kappa_diff(ii,jj)=(f_rho_3(jj)*10^(-8))/Zc^5;
                        kappa_grad(ii,jj)=kappa_diff(ii,jj)*(4.1868/10^(-2))+kappa_EU(ii,jj);
                 end

                elseif strcmp(visc,'Sutherland') 
        mu_grad(ii,jj) = mu_ref * (T_ref + S_mu_ref) ./ (T(ii) + S_mu_ref) .* (T(ii) / T_ref).^(3.0 / 2.0);
        kappa_grad(ii,jj) = kappa_ref * (T_ref + S_kappa_ref) ./ (T(ii) + S_kappa_ref) .* (T(ii) / T_ref).^(3.0 / 2.0);
                    
                elseif strcmp(visc,'constant')
		mu_grad(ii,jj)=0;
		kappa_grad(ii,jj)=0;
                end
      end
      end     
      % Rho
      for ii=1:length(T)
        mu_gradRho(ii)=  ddx_coeff*mu_grad(ii,:)'/(2*dRho);
        mu_gradRho2(ii)= d2dx2_coeff*mu_grad(ii,:)'/(dRho*dRho);
        kappa_gradRho(ii)=  ddx_coeff*kappa_grad(ii,:)'/(2*dRho);
        kappa_gradRho2(ii)= d2dx2_coeff*kappa_grad(ii,:)'/(dRho*dRho);
      end
      % T
      for jj=1:length(Rho)
        mu_gradT(jj)=  ddx_coeff*mu_grad(:,jj)/(2*dT);
        mu_gradT2(jj)= d2dx2_coeff*mu_grad(:,jj)/(dT*dT); 
        kappa_gradT(jj)=  ddx_coeff*kappa_grad(:,jj)/(2*dT);
        kappa_gradT2(jj)= d2dx2_coeff*kappa_grad(:,jj)/(dT*dT);
      end
      
      % Viscosity      
      mu_dim(i)=mu_grad(2,2);
      mu_T_dim(i)=mu_gradT(2);
      mu_TT_dim(i)=mu_gradT2(2);     
      mu_Rho_dim(i)=mu_gradRho(2);
      mu_RhoRho_dim(i)=mu_gradRho2(2);
      
      mu_RhoT_dim(i)=ddx_coeff*mu_gradRho(:)/(2*dT);
      mu_TRho_dim(i)=ddx_coeff*mu_gradT(:)/(2*dRho);
      
      % Thermal conductivity
      kappa_dim(i)=kappa_grad(2,2);
      kappa_T_dim(i)=kappa_gradT(2);
      kappa_TT_dim(i)=kappa_gradT2(2);
      kappa_Rho_dim(i)=kappa_gradRho(2);
      kappa_RhoRho_dim(i)=kappa_gradRho2(2);
      
      kappa_RhoT_dim(i)=ddx_coeff*kappa_gradRho(:)/(2*dT);
      kappa_TRho_dim(i)=ddx_coeff*kappa_gradT(:)/(2*dRho);
    
   end
   
   % Non-dimensional values with respect to infinity
   
  P_r=Zc/Rho_inf/T_inf/Ec_inf/Cp_inf;
  T_r=T_inf;
  Rho_r=Rho_inf;
  E_r=Zc/T_inf/Ec_inf/Cp_inf;
   
     % Pressure
  p=p_dim*P_r;
  p_T=p_T_dim*P_r*T_r; % T_inf
  p_TT=p_TT_dim*P_r*T_r^2;     
  p_Rho=p_Rho_dim*P_r*Rho_r;
  p_RhoRho=p_RhoRho_dim*P_r*Rho_r^2;

  p_RhoT=p_RhoT_dim*P_r*Rho_r*T_r;
  p_TRho=p_TRho_dim*P_r*Rho_r*T_r;

  % Viscosity    
  if strcmp(visc,'RG')
  mu=mu_dim/mu_inf;
  elseif strcmp(visc,'JST')
  mu=mu_dim/mu_inf;
  elseif strcmp(visc,'Sutherland')
  mu=mu_dim/mu_inf;
  elseif strcmp(visc,'constant')
  mu=mu0;  
  end
  
  mu_T=mu_T_dim/mu_inf*T_r; 
  mu_TT=mu_TT_dim/mu_inf*T_r^2; 
  mu_Rho=mu_Rho_dim/mu_inf*Rho_r; 
  mu_RhoRho=mu_RhoRho_dim/mu_inf*Rho_r^2;

  mu_RhoT=mu_RhoT_dim/mu_inf*Rho_r*T_r; 
  mu_TRho=mu_TRho_dim/mu_inf*Rho_r*T_r; 

  % Thermal conductivity
  if strcmp(visc,'RG')
  kappa=kappa_dim/kappa_inf;
  elseif strcmp(visc,'JST')
  kappa=kappa_dim/kappa_inf;
  elseif strcmp(visc,'Sutherland')
  kappa=kappa_dim/kappa_inf;
  elseif strcmp(visc,'constant')
  kappa=kappa0;
  end

  kappa_T=kappa_T_dim/kappa_inf*T_r; %%
  kappa_TT=kappa_TT_dim/kappa_inf*T_r^2; %%
  kappa_Rho=kappa_Rho_dim/kappa_inf*Rho_r; %%
  kappa_RhoRho=kappa_RhoRho_dim/kappa_inf*Rho_r^2;

  kappa_RhoT=kappa_RhoT_dim/kappa_inf*Rho_r*T_r; 
  kappa_TRho=kappa_TRho_dim/kappa_inf*Rho_r*T_r;

  % Internal energy
  e=e_dim*E_r;
  e_T=e_T_dim*E_r*T_r;
  e_Rho=e_Rho_dim*E_r*Rho_r;

  % Temperature
  PT_Rho=T_Rho_dim*T_r*Rho_r;
  RhoT_p=T_P_dim*T_r*P_r;

   if dT==dRho
    epsilon=dT;
   else 
    disp('dT is not equal to dRho!');
   end
end
