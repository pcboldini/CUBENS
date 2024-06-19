## Testcases
To successfully run CUBENS, copy the given input files to the corresponding folders and check output results

### Taylor-Green-Vortex (Sec. 5.1)
Run the 3-D TGV problem for:
  - `Re=1600, M=0.1, Pr=0.75`
  - `N=128^3, CFL=0.5`

input: `config_TGV.h` file with `Makefile`

output: `TGV_out.txt` file (time, Ekin, EInt, Enst)

### Boundary Layer  (Sec. 5.5.1)
Run the 2-D boundary layer for:
  - `Peng-Robinson EoS, Chung TP`
  - `Ec=0.01, M=0.203, Pr=2.395, T_red=0.9207, P_red=1.0844, Rho_red=2.2169`
  - no disturbance strip

input: `config_BL.h` file with `Makefile`, for inital condition: `main_PR.ipynb` file generating files in `inputDNS`

output: `ypl.1.` planes of non-conservative variables

### H-type breakdown (Sec. 5.2)
Run the initial stage of the H-type breakdown:
  - `Ideal gas EoS, Sutherland's law`
  - `M=0.2, Ec=0.016, Pr=0.75`
  - with 3-D disturbance strip

input: `config_BL.h` file with `Makefile`, for inital condition: `main_IG.ipynb` file generating files in `inputDNS`

output: `ruvwe.` restart file at `t=650000`


