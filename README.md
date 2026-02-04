<img src="CUBENS_logo.png" width="40%" height="40%">

## Main
CUBENS (CUBic Equation of state Navier-Stokes) is a massively-parallel GPU-accelerated high-order solver for direct numerical simulations of non-ideal wall-bounded flows. 
### Features
It incorporates:
 - non-ideal, strongly non-linear thermodynamics,
 - a wall-normal buoyant force,
 - high-order finite-difference schemes with convective terms in split and kinetic-energy- and entropy-preserving form,
 - pressure-equilibrium-preserving scheme,
 - non-ideal non-reflecting boundary conditions,
 - GPU-acceleration using OpenACC directives for computation offloading and asynchronous CUDA-aware MPI for GPU-GPU communication.

The following geometries can be simulated:
 - transitional boundary layer,
 - turbulent boundary layer,
 - channel flow.

Current developments:
 - suitable shock-capturing techniques for the supersonic and hypersonic flow regime,
 - curvilinear coordinates for more realistic geometries,
 - developing a new recycling/rescaling method for boundary layers with strongly non-ideal fluids.

## Motivation
The application of non-ideal fluids has rapidly and widely increased the relevance of a new branch of fluid mechanics called non-ideal compressible fluid dynamics. With respect to the rising number of industrial applications operating at non-ideal gas conditions, such as turbomachinery and heat exchangers, the development of more accurate theoretical, experimental, and numerical tools is required. In particular, one major challenge is the lack of knowledge in transitional and turbulent boundary layers due to the difficulty in performing experiments at high density and temperature conditions. On the contrary, high-fidelity simulations can significantly contribute to improve and accelerate the design of new engineering systems operating at non-ideal conditions. Fluids above their critical point play a key role in future energy conversion systems.

## Contributing

Any contributions and feedback that can improve CUBENS are appreciated. If you wish to contribute to the tool, please get in touch with the maintainers or open an Issue in the repository / a thread in Discussions.

### Reference 
For more information on CUBENS:

P.C. Boldini, R. Hirai, P. Costa, J.W.R. Peeters, R. Pecnik, "CUBENS: a GPU-accelerated high-order solver for wall-bounded flows with non-ideal fluids", Comput. Phys. Commun. (2025) 109507. [link_to_paper](https://www.sciencedirect.com/science/article/pii/S0010465525000104)

P.C. Boldini, "CUBENS: A GPU-accelerated high-order solver for wall-bounded flows with non-ideal fluids", Computer Physics Communications Seminar Series, [link_to_seminar](https://doi.org/10.52843/cassyni.gwmx6m)

## News
## **[2026/02/04]** Turbulent channel flow testcase added.
## **[2025/08/29]** Redlich-Kwong EoS is implemented.
# **[2025/01/17]** CUBENS v1.0 is finally online!
## **[2025/01/17]** The paper has been published (see link)!
## **[2025/01/14]** The paper has been accepted!

# Usage
## Compiling

1. CUBENS requires one of the following compilers :
 - for CPU architectures: GNU, Cray, Intel
 - for GPU architectures: NVHPC (NVIDIA), Cray (AMD)

 In the Makefile, the flag `ARCH` needs to be set

2. CUBENS requires the geometry input:
 - Boundary layer (BL), Channel (CHA), Taylor-Green-Vortex (TGV)

  In the Makefile, the flag `CASE` needs to be set   

3. CUBENS requires the Equation of State input:
 - Ideal gas (IG), Van der Waals (VdW), Peng-Robinson (PR)

 In the Makefile, the flag `EOS_LAW` needs to be set   

4. CUBENS requires the Transport Properties input:
 - Power Law, Sutherland, JossiStielThodos, Chung

 In the Makefile, the flag `VISC_LAW` needs to be set  

5. In the Makefile, the flag `BENCH` can be set for benchmark mode without parameter output
   
6. In the Makefile, the flag `DECOMP_OPTIONS` is set to double precision by default

7. If FFT wants to be performed, flag `FFT` needs to be active with the respective module installed on the system

To compile:
```
make
```
ATTENTION: each cluster needs specific modules! Please take a look into the cluster-specific requirements in case of errors.

## Running
The `config.h` file sets the simulation parameters

CUBENS can be easily run with
```
mpirun -np n_procs ./simulate n_rows n_cols 
```
where `n_procs` is the number of processors, `n_rows` is the number of processors in rows (j-direction), and `n_cols` is the number of processors in colums (k-direction)

NOTE: for `CASE = BL` and `CHA`, the initial condition in folder `preproc` has to be compiled a priori
 1. initBL: calculation of the compressible laminar solution and flow parameters -> `inputDNS/...`
 2. initBL: calculation of the flow parameters -> `inputDNS/...`

NOTE: for recycling-rescaling method `BC_inl_rescale = true`, the mean profiles have to be present in `preproc/turbRR`

On GPUs, the distribution of MPI and accelerators depends on the cluster architecture. An example for 8 GPUs on the Dutch National Supercomputer (Snellius, SURF), where one node is composed of 72 CPUs and 4 GPUs:
```
# Define partitions
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --gpus-per-node=4
#SBATCH --partition=gpu

# Load modules for MPI and other parallel libraries
module purge
module load 2022
module load foss/2022a NVHPC/22.7

# Execute the program in parallel 
mpirun -np 8 ./simulate 2 8 > out_sim
```

## Post-processing
The post-processing in CUBENS can be compiled with
```
make post
```
and can be run with
```
mpirun -np n_procs ./simulate 1 n_cols 
```
where the number of n_cols is fixed at 1. Note that the MPI distribution has to match with the `config.h` file. The results of post-processing are stored in go to `postproc/results`

For 2-D and 3-D plots, go to `postproc/visualize` (.xmf files).

## Interpolation
The interpolation for e.g. grid refinement in CUBENS can be compiled with
```
make interpol
```
and can be run with
```
mpirun -np n_procs ./interpol 1 n_cols 
```
where the number of n_cols is fixed at 1. Note that the MPI distribution has to match with the `config.h` file. New interpolated restart and planes are generated.

## Tests
See folder `/tests`.

## Final notes 

CUBENS is licensed under the MIT license.

