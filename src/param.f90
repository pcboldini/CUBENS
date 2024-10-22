! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! parameters module for the definition of all parameters of CUBENS

module mod_param
use decomp_2d
use mod_math
use io_std_units
implicit none 


! general parameters  
  real(mytype) :: Re      = 0.0_mytype
  real(mytype) :: delta99 = 0.0_mytype
  real(mytype) :: Pra     = 0.0_mytype
  real(mytype) :: Ma      = 0.0_mytype
  real(mytype) :: Ec      = 0.0_mytype
  real(mytype) :: Ri_wall = 0.0_mytype 
  real(mytype) :: Ri_unit = 0.0_mytype 
  real(mytype) :: dpdz    = 0.0_mytype
  real(mytype) :: Tinf    = 0.0_mytype
  real(mytype) :: Pinf    = 0.0_mytype
  real(mytype) :: Rhoinf  = 0.0_mytype
  real(mytype) :: Uinf    = 0.0_mytype
  real(mytype) :: Cpinf   = 0.0_mytype
  real(mytype) :: SOSinf  = 0.0_mytype
  real(mytype) :: Tref    = 0.0_mytype
  real(mytype) :: Pref    = 0.0_mytype
  real(mytype) :: Rhoref  = 0.0_mytype
  real(mytype) :: Uref    = 0.0_mytype
  real(mytype) :: Rref    = 0.0_mytype
  real(mytype) :: Cpref   = 0.0_mytype
  real(mytype) :: SOSref  = 0.0_mytype


! equation of state
#if defined(BL) || defined(CHA)
  character(len=20) :: USE_EOS 
#endif
  real(mytype) :: eos_dof = 9.0_mytype
  real(mytype) :: eos_ac = 0.22_mytype
  real(mytype) :: eos_Rgas = 200.0_mytype
  real(mytype) :: eos_Ru  = 8.31451_mytype
  !$acc declare create(eos_Rgas)
  real(mytype) :: Tcrit = 304.128200_mytype
  real(mytype) :: Pcrit = 7377300.0_mytype
  real(mytype) :: Vcrit = 2.138580e-03_mytype
  ! ideal gas 
  real(mytype) :: ig_gam = 1.4_mytype
  !$acc declare create(ig_gam)
  ! Van der Waals 
  real(mytype) :: vdw_a   = 3.0_mytype
  real(mytype) :: vdw_b   = 1.0_mytype/3.0_mytype  
  real(mytype) :: vdw_Zc  = 3.0_mytype/8.0_mytype  
  ! Peng-Robinson
  real(mytype) :: pr_a   = 0.45724_mytype  
  real(mytype) :: pr_b   = 0.07780_mytype  
  real(mytype) :: pr_Zc  = 0.3112_mytype 
  ! tables (REFPROP)
  real(mytype) :: rp_Zc  = 0.274586376_mytype


! Transport properties
#if defined(BL) || defined(CHA)
 character(len=20) :: USE_VISC
#endif 
  ! Sutherland (ideal gas)
  real(mytype) :: Smuref = 111.0_mytype
  ! JST & Chung
  real(mytype) :: Muref   = 0.0_mytype
  real(mytype) :: Kref  = 0.0_mytype
  ! boundary conditions
  character(len=20) :: wall_bc
  character(len=30) :: BC_bot 
  real(mytype)      :: Twall_bot = 0.0_mytype
  real(mytype)      :: Twall_top = 0.0_mytype
  character(len=30) :: BC_top
  character(len=30) :: BC_inl 
  character(len=30) :: BC_out
  character(len=30) :: BC_span
  logical, dimension(3) :: perBC = (/.false.,.true.,.false./)
  complex(mytype) :: alphaLST=(0.061251_mytype,0.013715_mytype)
  real(mytype) :: epsilonLST=1.0_mytype
  real(mytype) :: sigm_outlet = 0.0_mytype
  real(mytype) :: flag_supersonic_outlet = 0.0_mytype
  logical      :: BC_inl_rescale = .false.
  real(mytype) :: z_recycle = 15.0_mytype
  real(mytype) :: delta_inl = 1.0_mytype
  ! perturbation-related parameters
  integer      :: pert_calc = 0
  real(mytype) :: omega1 = 0.0_mytype
  real(mytype) :: omega2 = 0.0_mytype
  real(mytype) :: beta0 = 1.0_mytype
  real(mytype) :: pert_zLen  = 50.0_mytype
  real(mytype) :: pert_Reend = 100.0_mytype
  real(mytype) :: pert_Restart=100.0_mytype
  real(mytype) :: pert_zMid  = 50.0_mytype
  real(mytype) :: pert_ReMid = 1000.0_mytype
  real(mytype) :: pert_zSig  = 7.0_mytype
  real(mytype) :: pert_ySig  = 7.0_mytype
  real(mytype), dimension(:), allocatable :: pert_ampl 
  real(mytype), dimension(:), allocatable :: pert_F 
  real(mytype), dimension(:), allocatable :: pert_beta 
  namelist /paramPerturbation/ pert_calc, pert_ampl, pert_F, pert_beta, &
                              pert_zLen, pert_ReMid , pert_zSig, pert_ySig


! timestepping
  real(mytype) :: dtMax            = 1.2e-1_mytype 
  real(mytype) :: CFL              = 0.8_mytype
  integer      :: nsteps           = 1000000
  integer      :: readRestartFile  = -10000
  integer      :: intvSaveRestart  = 10000
  integer      :: intvSavePlanes   = 100 
  integer      :: intvSaveStats    = 100
  integer      :: savePlanesAfter  = 1
  integer      :: saveRestartAfter = 1
  integer      :: saveStatsAfter   = 1
  integer      :: intvCalcCFL      = 100
  integer      :: intvPrint        = 10
  integer      :: intvReadParam    = 100
  logical      :: abortSimulation  = .false.
  namelist /paramTimestepping/ dtMax, CFL, nsteps, readRestartFile, intvSaveRestart, &
                         intvSavePlanes, intvSaveStats, savePlanesAfter, saveRestartAfter, saveStatsAfter, &
                         intvCalcCFL, intvPrint, intvReadParam, abortSimulation


! FFT values: when pert_calc=1    
  integer :: fft_samples = 50 
  integer :: fft_step    = 400 
  namelist /paramFFT/ fft_samples, fft_step
  ! computational grid 
  integer :: imax,jmax,kmax
  integer :: p_row = 1, &
             p_col = 1
  real(mytype) :: len_x = 20.0_mytype         ! wall-normal
  real(mytype) :: len_y = 1.0_mytype          ! spanwise
  real(mytype) :: Redelta_start = 1.0_mytype
  real(mytype) :: Redelta_end = 1.0_mytype
  real(mytype) :: zStartDNS = 0.0_mytype
  real(mytype) :: zEndDNS = 0.0_mytype
  real(mytype) :: len_z = 0.0_mytype
  character(len=30) :: xmesh_type = "equid" ! equid, non_equid
  real(mytype) :: gridStretchX = 5.0_mytype   ! wall-normal clustering 
  real(mytype) :: ReTau = 0.0_mytype ! inlet friction Reynolds number
  character(len=30) :: zmesh_type = "equid" ! equid, non_equid
  real(mytype) :: z_1 = 0.3_mytype ! location of first bump in z-direction
  real(mytype) :: z_2 = 0.95_mytype ! location of second bump in z-direction
  real(mytype) :: bumpz1 = 0.2_mytype ! width of first bump in z-direction
  real(mytype) :: bumpz2 = 0.2_mytype ! width of second bump in z-direction
  real(mytype) :: zplus_min = 10.0_mytype ! zplus in turbulent region 
  real(mytype) :: zplus_max = 30.0_mytype ! zplus in laminar region
  real(mytype) :: zpert_1 = 0.3_mytype ! location of first bump in z-direction
  real(mytype) :: zpert_2 = 0.95_mytype ! location of second bump in z-direction
  real(mytype) :: bumpzpert1 = 0.2_mytype ! width of first bump in z-direction
  real(mytype) :: bumpzpert2 = 0.2_mytype ! width of second bump in z-direction
  real(mytype) :: zpluspert_min = 10.0_mytype ! zplus in turbulent region 
  ! output planes: index planes
  integer, dimension(:), allocatable :: yi_plane
  integer, dimension(:), allocatable :: xi_plane
  integer, dimension(:), allocatable :: zi_plane
  ! sponge: inlet, outlet, top
  real(mytype) :: spInlLen = 0.0_mytype
  real(mytype) :: spInlStr = 0.5_mytype
  real(mytype) :: spInlExp = 2.0_mytype
  real(mytype) :: spInLen_Reend = 100.0_mytype
  real(mytype) :: spOutLen = 0.0_mytype
  real(mytype) :: spOutStr = 0.5_mytype
  real(mytype) :: spOutExp = 2.0_mytype
  real(mytype) :: spOutLen_Resta = 100.0_mytype
  real(mytype) :: spTopLen = 0.0_mytype
  real(mytype) :: spTopStr = 0.5_mytype
  real(mytype) :: spTopExp = 2.0_mytype
  namelist /paramSponge/ spInlLen, spInlStr, spInlExp, &
                         spOutLen, spOutStr, spOutExp, &
                         spTopLen, spTopStr, spTopExp


! numerics
  integer      :: nStencilConv = 3
  integer      :: nStencilVisc = 2
  integer      :: nHalo = 1
  character(len=30) :: keep_flag = "classic" ! classic, pep


! postprocessing (postpro)
  integer :: p_row_pp  = 6 
  integer :: p_col_pp  = 72 
  character(len=30) :: avg_flag
  integer :: istart_pp = 147000
  integer :: iend_pp   = 246500
  integer :: istep_pp  = 500
  integer :: dt_step   = 0.001
  integer, dimension(:), allocatable :: stats_step             
  real(mytype), dimension(:), allocatable :: stats_time_rate
  integer :: rms_flag   = 1                                    
  integer :: istart_rms = 147000                               
  integer :: iend_rms   = 246500                               
  integer :: istep_rms  = 500                                  
  integer, dimension(:), allocatable :: index_rms_xpl          
  integer, dimension(:), allocatable :: index_rms_zpl          

  integer :: fft_flag  = 1
  integer, dimension(:), allocatable :: index_fft_span
  integer :: istart_fft = 147000                            
  integer :: iend_fft   = 246500                               
  integer :: istep_fft  = 500                                 
  integer :: index_fft_xpl = 15                             


! interpolation (interpol)
  integer :: timeStepRead = 0 ! step reading
  integer :: timeStepSave = 10 ! step writing
  integer :: nHaloInterpol = 0 ! new number of halo points
  integer :: inew = 100, & ! new number of points in x-direction
             jnew = 60, & ! new number of points in y-direction
             knew = 440 ! new number of points in z-direction
  character(len=30) :: xmesh_type_new = "equid" ! equid, non_equid
  real(mytype) :: gridStretchX_new = 4.0_mytype ! new grid stretching in x
  real(mytype) :: ReTau_new = 100.0_mytype
  character(len=30) :: zmesh_type_new = "equid" ! equid, non_equid
  real(mytype) :: z_1_new = 0.3_mytype ! location of first bump in z-direction
  real(mytype) :: z_2_new = 0.95_mytype ! location of second bump in z-direction
  real(mytype) :: bumpz1_new = 0.2_mytype ! width of first bump in z-direction
  real(mytype) :: bumpz2_new = 0.2_mytype ! width of second bump in z-direction
  real(mytype) :: zplus_min_new = 10.0_mytype ! zplus in turbulent region 
  real(mytype) :: zplus_max_new = 30.0_mytype ! zplus in laminar region
  real(mytype) :: len_x_new = 20.0_mytype, &  ! new length in x-direction
                  len_y_new = 7.0_mytype, & ! new length in y-direction
                  len_z_new = 200.0_mytype ! new length in z-direction
  real(mytype) :: Redelta_end_new = 1.0_mytype ! new end local Reynolds number in x-direction
  real(mytype) :: zEndDNS_new = 0.0_mytype ! new end domain in z-direction               
  integer :: p_row_intp = 1, &
             p_col_intp = 1
  integer, dimension(:), allocatable :: yi_plane_new
  integer, dimension(:), allocatable :: xi_plane_new
  integer, dimension(:), allocatable :: zi_plane_new
contains  
   

! this subroutine reads the config.h file and overwrites the simulation variables                               
  subroutine read_config()
  ! the number of halo cells is equal to the stencil of the convective fluxes
  nhalo=nStencilConv
#if defined(BL)
#include "../config_BL.h"
#elif defined(CHA)
#include "../config_CHA.h"
#elif defined(TGV)
#include "../config_TGV.h"
#endif
  end subroutine     


! this subroutine reads the init*_params.h file and overwrites the free-stream simulation variables                          
subroutine read_init_params()
#if defined(BL)
#include "../preproc/initBL/inputDNS/initBL_params.h"
  ! Streamwise parameters for boundary layer
  Re = Redelta_start * delta99 ! Reynolds_DNS
  zStartDNS = Redelta_start**2/Re ! zStart
  zEndDNS = Redelta_end**2/Re ! zEnd
  len_z = zEndDNS - zStartDNS ! streamwise length
  zEndDNS_new = Redelta_end_new**2/Re ! zEnd_new
  len_z_new = zEndDNS_new - zStartDNS ! new streamwise length
#elif defined(CHA)
#include "../preproc/initCHA/inputDNS/initCHA_params.h"
#endif
end subroutine


! write various parameters in params_variation.txt for real-time parameter variation
  subroutine io_writeParams_variation()
    use decomp_2d
    implicit none
    if (nrank == 0) then 
      open(17,file = 'postproc/params_variation.txt')
      write(17, nml=paramTimestepping)
#if defined(BL)
      write(17, nml=paramPerturbation)
      write(17, nml=paramSponge) 
#endif
      close(17)
    endif
  end subroutine


! read various parameters in params_variation.txt for real-time parameter variation
  subroutine io_readParams_variation()
    use mpi
    use decomp_2d
    implicit none
    integer :: ierr
    open(17,file = 'postproc/params_variation.txt')
    read(17, nml=paramTimestepping)
#if defined(BL)
    read(17, nml=paramPerturbation)
    read(17, nml=paramSponge) 
#endif
    close(17)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end subroutine


! print parameters on screen                                                                                                              
    subroutine print_init_params()
    use decomp_2d
    implicit none
    if (nrank == 0) then
      ! Dimensionless numbers
#if defined(BL)
      write(stdout,* ) 'BOUNDARY LAYER'
#elif defined(CHA)
      write(stdout,* ) 'CHANNEL'
#elif defined(TGV)
      write(stdout,* ) 'TAYLOR-GREEN VORTEX'
#elif defined(1D)
      write(stdout,* ) '1D TEST WAVE'
#endif
      write(stdout,* ) 
      write(stdout,* ) 'Dimensionless parameters:'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout, '(A, F10.4)') 'Mach number:                          ', Ma
      write(stdout, '(A, F10.4)') 'Reynolds number:                      ', Re
      write(stdout, '(A, F10.4)') 'Eckert number:                        ', Ec
      write(stdout, '(A, F10.4)') 'Prandtl number:                       ', Pra
      write(stdout, '(A, F10.4)') 'Pref:                                 ', Pref
      write(stdout, '(A, F10.4)') 'Wall Richardson number:               ', Ri_wall
      write(stdout, '(A, F10.4)') 'Richardson number (unit values):      ', Ri_unit
#if defined(BL)
      write(stdout, '(A, F10.4)') 'Inlet friction Reynolds number:       ', ReTau
#endif
#if defined(CHA)
      write(stdout, '(A, F10.4)') 'Initial streamwise pressure gradient: ', dpdz
#endif
#if defined(BL)
      write(stdout, '(A, F10.4)') 'Reynolds number end (delta):          ', Redelta_end
      write(stdout, '(A, F10.4)') 'delta_99/delta_blasius                ', delta99
#endif
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,* ) 
      ! EOS and FLUID
#if defined(IG)
      write(stdout,* ) 'EoS: ideal gas'
      write(stdout,* ) 
      write(stdout,* ) 'EoS parameters:'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,'(A, F10.4)') 'DOF:                                  ', eos_dof
      write(stdout,'(A, F10.4)') 'Gamma:                                ', ig_gam
      write(stdout,'(A, F10.4)') 'Rgas:                                 ', eos_Rgas
#elif defined(VdW)
      write(stdout,* ) 'EoS: Van der Waals'
      write(stdout,* ) 
      write(stdout,* ) 'EoS parameters:'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,'(A, F10.4)') 'DOF:                                  ',eos_dof
      write(stdout,'(A, F10.4)') 'Parameter a in VdW:                   ',vdw_a
      write(stdout,'(A, F10.4)') 'Parameter b in VdW:                   ',vdw_b
      write(stdout,'(A, F10.4)') 'Compressibility factor at crit. point:',vdw_Zc 

#elif defined(PR)
      write(stdout,* ) 'EoS: Peng-Robinson'
      write(stdout,* ) 
      write(stdout,* ) 'EoS parameters:'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,'(A, F10.4)') 'DOF:                                  ',eos_dof
      write(stdout,'(A, F10.4)') 'Parameter a in PR:                    ',pr_a
      write(stdout,'(A, F10.4)') 'Parameter b in PR:                    ',pr_b
      write(stdout,'(A, F10.4)') 'Compressibility factor at crit. point:',pr_Zc 
#endif
#if defined(VdW) || defined(PR)
      write(stdout,'(A, F10.4)') 'Acentric factor:                      ',eos_ac   
      write(stdout,'(A, F10.4)') 'Rgas [J/kg/K]:                        ',eos_Rgas 
      write(stdout,'(A, F10.4)') 'Tcrit [K]:                            ',Tcrit 
      write(stdout,'(A, F10.4)') 'Pcrit [bar]:                          ',Pcrit*1e-5 
      write(stdout,'(A, F10.4)') 'Vcrit [m3]:                           ',Vcrit 
      write(stdout,'(A, F10.4)') 'Reduced Temperature [-]:              ',Tref
      write(stdout,'(A, F10.4)') 'Reduced Pressure [-]:                 ',Pref
      write(stdout,'(A, F10.4)') 'Reduced Density [-]:                  ',Rhoref
      write(stdout,'(A, F10.4)') 'Specific gas constant [-]:            ',Rref
      write(stdout,'(A, F10.4)') 'Heat capacity cp [-]:                 ',Cpref
      write(stdout,'(A, F10.4)') 'Speed of sound [-]:                   ',SOSref 
#endif
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,* )  
      ! Viscosity and conductivity  
#if defined(PowerLaw)
      write(stdout,* ) 'Transport properties: Power Law' 
      write(stdout,* ) 
      write(stdout,* ) 'TP parameters:'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,'(A, F10.4)') 'Power exponent:                       ',3.0_mytype/4.0_mytype 
#elif defined(Sutherland)
      write(stdout,* ) 'Transport properties: Sutherland' 
      write(stdout,* ) 
      write(stdout,* ) 'TP parameters:'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,'(A, F10.4)') 'Sutherland constant [K]:              ',Smuref
      write(stdout,'(A, F10.4)') 'Reference temperature [K]:            ', Tinf 
#elif defined(JST)  
      write(stdout,* ) 'Transport properties: Jossi, Stiel and Thodos' 
      write(stdout,* ) 
      write(stdout,* ) 'TP parameters:'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,'(A, F10.4)') 'Ref. viscosity [-]:                   ',Muref
      write(stdout,'(A, F10.4)') 'Ref. conductivity [-]:                ',Kref 
#elif defined(Chung)  
      write(stdout,* ) 'Transport properties: Chung' 
      write(stdout,* ) 
      write(stdout,* ) 'TP parameters:'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,'(A, F10.4)') 'Ref. viscosity [Pas]:                 ',Muref
      write(stdout,'(A, F10.4)') 'Ref. conductivity [W/m/K]:            ',Kref 
#endif
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,* )
    endif
  end subroutine

  
end module mod_param
