function[xSol,ySol,CSol,KappaSol,PrSol,RhoSol,TSol_dim] = PR_RK4_BF(dEqs,y,PR_PROP,BL_PARA,N,dh)
 % Initial values   
x=0;
xSol = zeros(N,numel(x));
ySol = zeros(N,numel(y));
CSol = zeros(N,1);
PrSol = zeros(N,1);
RhoSol = zeros(N,1);
KappaSol = zeros(N,1);
TSol_dim = zeros(N,1);
k = 1;
h_ref=y(4);

% Evaluation of y at the wall
xSol(1) = x; ySol(1,:) = y;
% Evaluation of C1 and Pr at the wall
[CSol(1),KappaSol(1),PrSol(1),RhoSol(1),TSol_dim(1)]=Calc_PR_Prop(h_ref,PR_PROP,BL_PARA);
for i=1:N-1
        C=CSol(i); Pr=PrSol(i); Rho=RhoSol(i);
        k1=feval(dEqs,x,y,C,Pr,Rho);
        k2=feval(dEqs,x+dh(i)/2,y+k1*dh(i)/2,C,Pr,Rho);
        k3=feval(dEqs,x+dh(i)/2,y+k2*dh(i)/2,C,Pr,Rho);
        k4=feval(dEqs,x+dh(i),y+k3*dh(i),C,Pr,Rho);
        dy=(dh(i)/6)*(k1+2*k2+2*k3+k4);
        x = x+dh(i);
        y = y + dy;
        k = k+1;
        xSol(k) = x; ySol(k,:) = y;  
        [CSol(k),KappaSol(k),PrSol(k),RhoSol(k),TSol_dim(k)]=Calc_PR_Prop(ySol(k,4),PR_PROP,BL_PARA);        
end
end

