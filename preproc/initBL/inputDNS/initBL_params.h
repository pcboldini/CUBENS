!------------ use equation of state ------------------ 
USE_EOS  = "VdW" 
!------------ set non-dimensional free-stream values for computation ------------------ 
delta99   = 4.4979046926 
Pra   = 1.0000000000 
Ec   = 5.0000000000e-01 
Ma   = 1.1219083278e+00 
eos_dof   = 9.000000 
eos_ac   = 0.224000 
!------------ set dimensional free-stream values for computation ------------------ 
Tcrit   = 304.128200 
Pcrit   = 7377300.000000 
Vcrit   = 2.920638e-03 
eos_Rgas   = 188.924058 
! ------------ set wall BC for computation ------------------ 
wall_bc  = "adiab" 
Twall_bot   = -1.000000 
!------------ set reference free-stream values for computation ------------------ 
Tref   = 0.9000000000 
Pref   = 1.1000000000 
Rhoref   = 1.8047114354 
Cpref   = 8.0239456825 
SOSref   = 2.7658411983 
Rref   = 2.6666666667 
! ----------- set viscosity and conductivity ---------------- 
USE_VISC  = "JST" 
Muref   = 1.6400035224e-03 
Kref   = 3.5314649926e-02 
