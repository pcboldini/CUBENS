!------------ use equation of state ------------------ 
USE_EOS  = "VdW" 
!------------ set non-dimensional free-stream values for computation ------------------ 
delta99   = 5.0819683966 
Pra   = 1.0000000000 
Ec   = 3.1141000000e-03 
Ma   = 2.0000054386e-01 
eos_dof   = 9.000000 
eos_ac   = 0.224000 
!------------ set dimensional free-stream values for computation ------------------ 
Tcrit   = 304.128200 
Pcrit   = 7377300.000000 
Vcrit   = 2.920638e-03 
eos_Rgas   = 188.924058 
! ------------ set wall BC for computation ------------------ 
wall_bc  = "isoth" 
Twall_bot   = 1.090909 
!------------ set reference free-stream values for computation ------------------ 
Tref   = 1.1000000000 
Pref   = 1.1000000000 
Rhoref   = 0.5801061435 
Cpref   = 8.8869952124 
SOSref   = 1.4246011232 
Rref   = 2.6666666667 
! ----------- set viscosity and conductivity ---------------- 
USE_VISC  = "JST" 
Muref   = 5.4580934235e-04 
Kref   = 2.5270360528e-02 
