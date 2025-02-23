!------------ use equation of state ------------------ 
USE_EOS  = "IG" 
!------------ set non-dimensional free-stream values for computation ------------------ 
delta99   = 4.9982781155 
Pra   = 0.7200000000 
Ec   = 1.6000000000e-02 
Ma   = 2.0000000000e-01 
Pref   = 17.857143 
ig_gam   = 1.400000 
eos_Rgas   = 17.857143 
eos_dof   = 9.000000 
Ri_wall   = -0.010000 
Ri_unit   = 0.210000 
! ------------ set wall BC for computation ------------------ 
wall_bc  = "isoth" 
Twall_bot   = 1.050000 
! ----------- set viscosity and conductivity ---------------- 
USE_VISC  = "Sutherland" 
Smuref   = 111.000000 
Tinf   = 300.000000 
