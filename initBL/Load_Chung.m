function [PROP] = Load_Chung(PROP)

omega_ac=PROP.omega_ac;
dof=PROP.dof;

PROP.A=1.16145; PROP.B=0.14874; PROP.C=0.52487; PROP.D=0.77320; PROP.E=2.16178; PROP.F=2.43787; PROP.G=-6.435*10^(-4);
PROP.H=7.27371; PROP.S=18.0323; PROP.W=-0.76830;

PROP.Fc=1-0.2756*omega_ac; % only non-polar

PROP.a0=[6.32402,0.0012102,5.28346,6.62263,19.7454,-1.89992,24.2745,0.79716,-0.23816,0.068629];
PROP.a1=[50.4119,-0.0011536,254.209,38.0957,7.63034,-12.5367,3.44945,1.11764,0.067695,0.34793];
for i=1:numel(PROP.a0)
    PROP.Arho(i)=PROP.a0(i)+PROP.a1(i)*omega_ac;
end

PROP.b0=[2.41657,-0.50924,6.61069,14.54250,0.79274,-5.86340,81.171];
PROP.b1=[0.74824,-1.50936,5.62073,-8.91387,0.82019,12.80050,114.158];
for i=1:numel(PROP.b0)
    PROP.Brho(i)=PROP.b0(i)+PROP.b1(i)*omega_ac;
end

PROP.alpha=dof/2-3/2;
PROP.beta=0.7862-0.7109*omega_ac+1.3168*omega_ac^2; % only non-polar

end