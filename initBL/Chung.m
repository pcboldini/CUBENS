function [mu,kappa] = Chung(rho,T,PROP)
global T_crit M

% Loading 

v_crit=PROP.v_crit;
Zc_RP=PROP.Zc_RP;
Zc=PROP.Zc;
M=PROP.M;
T_crit=PROP.T_crit;
T_dim=T*T_crit;

% Constants 
A=PROP.A; B=PROP.B; C=PROP.C; D=PROP.D; E=PROP.E; F=PROP.F; G=PROP.G;
H=PROP.H; S=PROP.S; W=PROP.W;
Fc=PROP.Fc;
Arho=PROP.Arho;
Brho=PROP.Brho;
alpha=PROP.alpha;
beta=PROP.beta;

% Viscosity 

V_crit=v_crit*M*10^6;
ek=T_crit/1.2593;
T_dimless=T_dim./ek;
psi=A./T_dimless.^(B)+C./exp(D.*T_dimless)+E./exp(F.*T_dimless)+G.*T_dimless.^(B).*sin(S.*T_dimless.^W-H);

mu_0=4.0785 * 10^(-6)* Fc*( M* 10^3 * T_dim ).^0.5 ./ ( V_crit^(2/3) * psi ); % new version

Y=rho*Zc_RP/Zc/6;
G1=(1-0.5.*Y)./(1-Y).^3;
G2=( Arho(1).*(1-exp(-Arho(4).*Y))./Y + Arho(2).*G1.*exp(Arho(5).*Y) + Arho(3).*G1 )./( Arho(1).*Arho(4)+Arho(2)+Arho(3) );

mu_k=mu_0.*(1./G2+Arho(6).*Y);

mu_p=3.6344*10^(-6)* ( M* 10^3 * T_crit ).^0.5 ./ ( V_crit^(2/3) )...
    .*( Arho(7).*Y.^2.*G2.*exp(Arho(8)+Arho(9).*T_dimless.^(-1)+Arho(10).*T_dimless.^(-2)) );

mu=mu_k+mu_p;

% Conductivity 

Z=2.0+10.5*T.^2;
Psi=1+alpha*( 0.215+0.28288*alpha-1.061*beta+0.26665.*Z )./( 0.6366 + beta*Z+1.061*alpha*beta);

kappa_0=31.2*(mu_0./(M)).*Psi;

H2=( Brho(1).*(1-exp(-Brho(4).*Y))./Y + Brho(2).*G1.*exp(Brho(5).*Y) + Brho(3).*G1 )./( Brho(1).*Brho(4)+Brho(2)+Brho(3) );

kappa_k=kappa_0.*(1./H2+Brho(6).*Y);

kappa_p=3.586*10^(-3)* (  T_crit / (M) ).^0.5 ./ ( V_crit^(2/3) )...
    .*( Brho(7).*Y.^2.*H2.*T.^(0.5) );

kappa=kappa_k+kappa_p;
end

