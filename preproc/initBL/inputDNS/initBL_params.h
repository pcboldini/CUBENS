!------------ use equation of state ------------------ 
USE_EOS  = "IG" 
!------------ set non-dimensional free-stream values for computation ------------------ 
delta99   = 4.9181043571 
Pra   = 0.7500000000 
Ec   = 4.0000000000e-03 
Ma   = 1.0000000000e-01 
Pref   = 71.428571 
ig_gam   = 1.400000 
eos_Rgas   = 71.428571 
eos_dof   = 9.000000 
! ------------ set wall BC for computation ------------------ 
wall_bc  = "adiab" 
Twall_bot   = -1.000000 
! ----------- set viscosity and conductivity ---------------- 
USE_VISC  = "Sutherland" 
Smuref   = 111.000000 
Tinf   = 300.000000 
