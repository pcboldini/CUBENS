!------------ use equation of state ------------------ 
USE_EOS  = "VdW" 
!------------ set non-dimensional reference values for computation ------------------ 
Pra   = 1.0000000000 
Ec   = 4.0000000000e-03 
Ma   = 1.0034653142e-01 
eos_dof   = 9.000000 
eos_ac   = 0.224000 
!------------ set dimensional values for computation ------------------ 
Tcrit   = 304.128200 
Pcrit   = 7377300.000000 
Vcrit   = 2.920638e-03 
eos_Rgas   = 188.924058 
! ------------ set wall BC for computation ------------------ 
Twall_bot   = 1.000000 
Twall_top   = 1.000000 
!------------ set reference values for computation ------------------ 
Tref   = 0.9000000000 
Pref   = 1.1000000000 
Rhoref   = 1.8047114354 
Cpref   = 8.0239456825 
SOSref   = 2.7658411983 
Rref   = 2.6666666667 
! ----------- set viscosity and conductivity ---------------- 
USE_VISC  = "Constant" 
