function [Re_BL,Cf,St] = extract_selfsim(Dy_BL,FLOW,THERMO,deltaHw,BL_PARA)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Re_BL=linspace(1e2,1e3,200)';
Ec_inf=BL_PARA.Ec_inf;
Pr_inf=BL_PARA.Pr_inf;

for i=1:numel(Re_BL)
    U(:,i)=FLOW(:,2);
    mu(:,i)=THERMO(:,1);
    T(:,i)=FLOW(:,3);
    kappa(:,i)=THERMO(:,2);

    Uy(:,i)=Dy_BL*U(:,i);
    Ty(:,i)=Dy_BL*T(:,i);
end

for i=1:numel(Re_BL)
    Sy(:,i)=mu(:,i).*Uy(:,i);
    Qy(:,i)=kappa(:,i).*Ty(:,i);
    Cf(i)=2*Sy(1,i)./Re_BL(i);
    St(i)=Qy(1,i)./deltaHw(i)/Re_BL(i)/Pr_inf/Ec_inf;
end

end

