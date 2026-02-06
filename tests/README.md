## Testcases
To successfully run CUBENS, copy the given input files to the corresponding folders and check output results

### Taylor-Green-Vortex (Sec. 5.1)
Run the 3-D TGV problem for:
  - `Re=1600, M=0.1, Pr=0.75`
  - `N=128^3, CFL=0.5`

input: `config_TGV.h` file with `Makefile`

output: `tests/TGV/output.txt` file (time, Ekin, EInt, Enst)

### H-type breakdown (Sec. 5.2)
Run the initial stage of the H-type breakdown:
  - `Ideal gas EoS, Sutherland's law`
  - `M=0.2, Ec=0.016, Pr=0.75`
  - with 3-D disturbance strip

input: `config_BL.h` file with `Makefile`, for inital condition: `main_IG.ipynb` file generating files in `preproc/initBL/inputDNS`

output: `ypl.1.` planes of non-conservative variables at `t=650000` (periodic solution)

### Turbulent Boundary Layer (Sec. 5.3)
Run the recycling-rescaling method of a turbulent boundary layer
  - `Ideal gas EoS, power law`
  - `M=0.2, Ec=0.016, Pr=0.75`
  - `Recycling position: x_rcy=112.3`
  - `Inlet friction Reynolds number: Re_tau,inl = 112.9`

input: `config_BL.h` file with `Makefile`, for inital mean-flow condition: `mean_values.txt` file `preproc/turbRR/`

output: time- and spanwise-average planes of streamwise velocity and shear stress at `t=400000`, wall-normal properties in `wall_normal_prop.txt`

### Laminar Boundary Layer  (Sec. 5.5.1)
Run the 2-D boundary layer for:
  - `Peng-Robinson EoS, Chung TP`
  - `Ec=0.01, M=0.203, Pr=2.395, T_r=0.9207, P_r=1.0844, Rho_r=2.2169`
  - no disturbance strip

input: `config_BL.h` file with `Makefile`, for inital condition: `main_PR.ipynb` file generating files in `preproc/initBL/inputDNS`

output: `ypl.1.` planes of non-conservative variables at `t=100000` (periodic solution)

### Turbulent Channel Flow (Extra)
Run the 3-D channel flow at `Re_tau=180` (Kim et al., JFM, 1987) for:
  - `Re=2800, M=0.2, Pr=0.75, dpdx=2.0`
  - `N_x=129, N_y=160, N_z=192, L_x=2.0, L_y=2.0\pi, L_z=4.0\pi`

input: `config_CHA.h` file with `Makefile`, for inital condition: `main_IG.ipynb` file generating files in `preproc/initCHA/inputDNS`

output: time- and spanwise-average profiles at `t=1200000` in `stats`, postprocessing script in `postproc.ipynb`, reference files in `literature`

### Differentially Heated Cavity (Extra)
Run the DHC problem at `Ra=10^6` (Le Quéré et al., ESAIM, 2005) for:
  - `Re=1532.13, M=0.2, Pr=0.71, 1/Fr^2=0.5`
  - `T_h=1.6, T_c=0.4, L=1.0, N_x=N_z=301, VISC_LAW=Constant`

input: `config_DHC.h` file with `Makefile`

output: `ypl.1.` planes and restart file of non-conservative variables at `t=2000000`, postprocessing script in `postproc.ipynb`