
!----------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                                  !
!      mod_param                                                                                                                   !
!                                                                                                                                  !
!      DESCRIPTION:                                                                                                                !
!      ------------                                                                                                                !
!      This module is for the input and output of the general configuration.                                                       !
!      The first part is just to initialise the variables                                                                          !
!                                                                                                                                  !
!                                                                                                                                  !
!      CHANGELOG:                                                                                                                  !
!      ----------                                                                                                                  !
!         xx.xx.2021: module created (Rene)                                                                                        !
!         03.12.2021: output for simulation parameters added (Pietro)                                                              !
!         07.12.2021: input free-stream paramters and viscosity added (Pietro)                                                     ! 
!                                                                                                                                  !
!----------------------------------------------------------------------------------------------------------------------------------!


module mod_param

use decomp_2d
use mod_math
use io_std_units

implicit none 

! ------------ initialization -------------------------------
! 
  character(len=40) :: CASE 
!
! ------------ set values for computation ------------------------
!  
  real(mytype) :: Re      = 0.0_mytype
  real(mytype) :: delta99 = 0.0_mytype
  real(mytype) :: Pra     = 0.0_mytype
  real(mytype) :: Ma      = 0.0_mytype
  real(mytype) :: Ec      = 0.0_mytype
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
!
! ------------ shock capturing ------------------------------
! 
  logical :: useArtificialBulkViscosity = .false. 
!
! ------------ Equation of state ----------------------------
!
  character(len=20) :: USE_EOS 
  real(mytype) :: eos_dof = 9.0_mytype
  real(mytype) :: eos_ac = 0.22_mytype
  real(mytype) :: eos_Rgas = 200.0_mytype
  real(mytype) :: eos_Ru  = 8.31451_mytype ! exact!
  !$acc declare create(eos_Rgas)

  real(mytype) :: Tcrit = 304.128200_mytype
  real(mytype) :: Pcrit = 7377300.0_mytype
  real(mytype) :: Vcrit = 2.138580e-03_mytype
!
! ------------ ideal gas ------------------------------------
!
  real(mytype) :: ig_gam = 1.4_mytype
  !$acc declare create(ig_gam)
!
! ------------ van der Waals --------------------------------
!
  real(mytype) :: vdw_a   = 3.0_mytype  !exact!
  real(mytype) :: vdw_b   = 1.0_mytype/3.0_mytype  !exact!
  real(mytype) :: vdw_Zc  = 3.0_mytype/8.0_mytype  !exact!
!
! ------------ Peng-Robinson --------------------------------
!
  real(mytype) :: pr_a   = 0.45724_mytype  !exact!
  real(mytype) :: pr_b   = 0.07780_mytype  !exact!
  real(mytype) :: pr_Zc  = 0.3112_mytype  !exact!
!
! ------------ REFPROP --------------------------------
!
  real(mytype) :: rp_Zc  = 0.274586376_mytype  !exact!
!
! ------------ visc and cond --------------------------------
!
 character(len=20) :: USE_VISC
!
! ------------ ideal gas ------------------------------------
!
 real(mytype) :: Stref  = 0.0_mytype
 real(mytype) :: Muinf  = 0.0_mytype
 real(mytype) :: Muref  = 0.0_mytype
 real(mytype) :: Smuref = 0.0_mytype
 real(mytype) :: Kinf   = 0.0_mytype

! ------------ JST ----------------------------------------
!
 real(mytype) :: Kref   = 0.0_mytype
 real(mytype) :: Skref  = 0.0_mytype
!
! ------------ Chung ----------------------------------------
!
! ------------ boundary conditions ------------------------------
!
!  adiab_std,  isoth_std  : standard implementation of wall BC
!  adiab_nrbc, isoth_nrbc : non-reflecting implementation of wall BC
!  freestream             : non-reflecting free stream BC
!  inlet_nrbc             : non-reflecting inlet BC, setting: rho,u,v,w 
!  inlet_std_subsonic     : standard implementation of inlet BC, setting: rho,u,v,w, extrapolating p 
!  inlet_std_supersonic   : standard implementation of inlet BC, setting: rho,u,v,w,p 
!
  character(len=20) :: wall_bc
  character(len=30) :: BC_bot = "isothermal_collocated"
  real(mytype)      :: Twall_bottom = 0.0_mytype
  real(mytype)      :: Twall_top = 0.0_mytype
  character(len=30) :: BC_top = "isothermal_collocated"
  character(len=30) :: BC_inl = "inlet_std_subsonic"
  character(len=30) :: BC_out = "outlet_std_subsonic"
  logical, dimension(3) :: perBC = (/.false.,.true.,.false./)

  complex(mytype) :: alphaLST=(0.061251_mytype,0.013715_mytype)
  real(mytype) :: epsilonLST=1.0_mytype
  real(mytype) :: sigm_outlet = 0.0_mytype
  real(mytype) :: flag_supersonic_outlet = 0.0_mytype
!
! ------------ Perturbation Related Parameters ------------------------
!
  ! integer      :: kC = 40   ! Central node of peturbation
  ! integer      :: LP = 0   ! Length of peturbation
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
!
! ------------ Wall heating Related Parameters ------------------------
!
  real(mytype) :: Twall_new = 100.0_mytype
  real(mytype) :: ReTem1 = 1.0_mytype
  real(mytype) :: ReTem2 = 1.0_mytype
  real(mytype) :: deltatem1 = 1.0_mytype
  real(mytype) :: deltatem2 = 1.0_mytype
!
! ------------ timestepping ------------------
!  
  real(mytype) :: dtMax           = 1.2e-1_mytype  ! 2*pi/omega/10000
  real(mytype) :: CFL             = 0.8_mytype
  integer      :: nsteps          = 1000000
  integer      :: readRestartFile = -10000
  integer      :: intvSaveRestart = 10000
  integer      :: intvSavePlanes  = 100 
  integer      :: savePlanesAfter = 1
  integer      :: saveRestartAfter = 1
  integer      :: intvCalcCFL     = 100
  integer      :: intvPrint       = 10
  integer      :: intvReadParam   = 100
  logical      :: abortSimulation = .false.

  namelist /paramTimestepping/ dtMax, CFL, nsteps, readRestartFile, intvSaveRestart, &
                         intvSavePlanes, savePlanesAfter, saveRestartAfter, intvCalcCFL, &
                         intvPrint, intvReadParam, abortSimulation
!
! -------------------------- values for FFT ------------------------         
!
  integer :: fft_samples = 50 ! per period
  integer :: fft_step    = 400 ! per sample

  namelist /paramFFT/ fft_samples, fft_step
!
! ------------ set values for computational grid  ------------------
!
  integer :: imax,jmax,kmax

  integer :: p_row = 1, &
             p_col = 1 !140 !24 !144

  integer :: nHalo = 3 ! currently only 3 is supported

  real(mytype) :: len_x = 20.0_mytype         ! wall-normal
  real(mytype) :: len_y = 1.0_mytype          ! spanwise
  real(mytype) :: Redelta_start = 1.0_mytype
  real(mytype) :: Redelta_end = 1.0_mytype
  real(mytype) :: zStartDNS = 0.0_mytype
  real(mytype) :: zEndDNS = 0.0_mytype
  real(mytype) :: len_z = 0.0_mytype

  character(len=30) :: xmesh_type = "equid" ! equid, non_equid
  real(mytype) :: gridStretchX = 5.0_mytype   ! clustering 
  real(mytype) :: ReTau = 0.0_mytype

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
!
! ------------ output planes ------------------
!    
  integer ::  yi_plane = 1
  integer ::  xi_plane = 25
  integer ::  zi_plane = 100
!
! --------------------- set values for sponges  ----------------------
!
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
!
! ------------------------------------------------------------------------------------------
! 
! ------ This are the parameters for postpro ------------------
!
  integer :: p_row_pp  = 6 
  integer :: p_col_pp  = 72 
  integer :: istart_pp = 147000
  integer :: iend_pp   = 246500
  integer :: istep_pp  = 500

  integer :: fft_flag  = 1
  integer :: rms_flag  = 1
  integer, dimension(:), allocatable :: index_fft_span


!
! parameters for interpolating a solution onto a new mesh 
!
  integer :: timeStepRead = 0
  integer :: timeStepSave = 10

  integer :: inew = 100, &
             jnew = 60, &
             knew = 440

  integer :: nHaloInterpol = 0

  character(len=30) :: xmesh_type_new = "equid" ! equid, non_equid
  real(mytype) :: gridStretchX_new = 4.0_mytype
  real(mytype) :: ReTau_new = 100.0_mytype

  character(len=30) :: zmesh_type_new = "equid" ! equid, non_equid
  real(mytype) :: z_1_new = 0.3_mytype ! location of first bump in z-direction
  real(mytype) :: z_2_new = 0.95_mytype ! location of second bump in z-direction
  real(mytype) :: bumpz1_new = 0.2_mytype ! width of first bump in z-direction
  real(mytype) :: bumpz2_new = 0.2_mytype ! width of second bump in z-direction
  real(mytype) :: zplus_min_new = 10.0_mytype ! zplus in turbulent region 
  real(mytype) :: zplus_max_new = 30.0_mytype ! zplus in laminar region
  
  real(mytype) :: len_x_new = 20.0_mytype, &
                  len_y_new = 7.0_mytype, &
                  len_z_new = 200.0_mytype  !100.0  

  real(mytype) :: Redelta_end_new = 1.0_mytype  
  real(mytype) :: zEndDNS_new = 0.0_mytype                      

  integer :: p_row_intp = 1, &
             p_col_intp = 1

  integer ::  yi_plane_new = 1
  integer ::  xi_plane_new = 25
  integer ::  zi_plane_new = 100

contains

!-------------------------------------------------------------------------------------------------------------------------------!
!   read_config                                                                                                                 !
!                                                                                                                               !
!   This subroutine reads the config.h file and overwrites the simulation variables                                             !
!                                                                                                                               !
!   CHANGELOG:                                                                                                                  !
!   ----------                                                                                                                  !
!      xx.xx.2021: subroutine created (Rene)                                                                                    !
!      02.12.2021: config file modified (Pietro)                                                                                !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!

  subroutine read_config()

#include "../config.h"

  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------!
!   read_initBL_params                                                                                                          !
!                                                                                                                               !
!   This subroutine reads the initBL_params.h file and overwrites the free-stream simulation variables                          !
!                                                                                                                               !
!   CHANGELOG:                                                                                                                  !
!   ----------                                                                                                                  !
!      xx.xx.2021: subroutine created (Rene)                                                                                    !
!      02.12.2021: initBL_params generated (Pietro)                                                                             !
!      09.02.2022: VDW added (Pietro)                                                                                           !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!

subroutine read_initBL_params()

#include "../initBL/inputDNS/initBL_params.h"

  Re = Redelta_start * delta99 ! Reynolds_DNS
  zStartDNS = Redelta_start**2/Re ! zStart
  zEndDNS = Redelta_end**2/Re ! zEnd
  len_z = zEndDNS - zStartDNS ! streamwise length

  zEndDNS_new = Redelta_end_new**2/Re ! zEnd_new
  len_z_new = zEndDNS_new - zStartDNS ! new streamwise length

end subroutine

  subroutine io_writeParams_variation()
    use decomp_2d
    implicit none
    if (nrank == 0) then 
      open(17,file = 'postproc/params_variation.txt')
      write(17, nml=paramTimestepping)
      write(17, nml=paramPerturbation)
      write(17, nml=paramSponge) 
      close(17)
    endif
  end subroutine


  subroutine io_readParams_variation()
    use mpi
    use decomp_2d
    implicit none
    integer :: ierr
    open(17,file = 'postproc/params_variation.txt')
    read(17, nml=paramTimestepping)
    read(17, nml=paramPerturbation)
    read(17, nml=paramSponge) 
    close(17)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------!
!   print_params                                                                                                                 !
!                                                                                                                               !
!   This subroutine prints all simulation parameters.                                                                           !
!                                                                                                                               !
!   CHANGELOG:                                                                                                                  !
!   ----------                                                                                                                  !
!      03.12.2021: subroutine created (Pietro)                                                                                  !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!

    subroutine print_init_params()
    use decomp_2d
    implicit none

    if (nrank == 0) then

      ! Dimensionless numbers
      write(stdout,* ) 
      write(stdout,* ) '-------------------------'
      write(stdout,* ) 'Flow case                                ',CASE
      write(stdout,* ) 
      write(stdout,* ) 'Dimensionless parameters:'
      write(stdout,* ) '-------------------------'
      write(stdout,* ) 'Mach number                              ',Ma
      write(stdout,* ) 'Reynolds number start (delta_99)         ',Re
      write(stdout,* ) 'Reynolds number end (delta)              ',Redelta_end
      write(stdout,* ) 'delta_99/delta_blasius                   ',delta99
      write(stdout,* ) 'Eckert number                            ',Ec
      write(stdout,* ) 'Prandtl number                           ',Pra
      write(stdout,* ) 'Pref                                     ',Pref
      write(stdout,* ) 'Reynolds tau (laminar inlet)             ',ReTau
      write(stdout,* ) '-------------------------'
      write(stdout,* ) 

      
      ! EOS and FLUID
      write(stdout,* ) 'EOS                                       ',USE_EOS
      write(stdout,* ) 
      write(stdout,* ) 'EoS parameters:'
      write(stdout,* ) '-------------------------'
      write(stdout,* ) 'DOF                                     ',eos_dof
      select case (USE_EOS)
      case("IG")
        write(stdout,* ) 'Gamma                                    ',ig_gam
        write(stdout,* ) 'Rgas                                     ',eos_Rgas
      case("VdW")   ! Van der Waals gas 
        write(stdout,* ) 'Parameter a in VdW                       ',vdw_a
        write(stdout,* ) 'Parameter b in VdW                       ',vdw_b
        write(stdout,* ) 'Compressibility factor at critical point ',vdw_Zc 
        write(stdout,* ) 'Acentric factor                          ',eos_ac       
      case("PR")   ! Peng-Robinson gas 
        write(stdout,* ) 'Parameter a in PR                        ',pr_a
        write(stdout,* ) 'Parameter b in PR                        ',pr_b
        write(stdout,* ) 'Compressibility factor at critical point ',pr_Zc 
        write(stdout,* ) 'Acentric factor                          ',eos_ac 
      end select
      select case (USE_EOS)
      case("IG") 
        write(stdout,* ) '-------------------------'
        write(stdout,* ) 
        write(stdout,* ) 'Dimensional parameters:'
        write(stdout,* ) '-------------------------'
        write(stdout,* ) 'Tinf                          [K]        ',Tinf 
      case("VdW","PR")   
        write(stdout,* ) '-------------------------'
        write(stdout,* ) 
        write(stdout,* ) 'Dimensional parameters:'
        write(stdout,* ) '-------------------------'
        write(stdout,* ) 'Rgas                                     ',eos_Rgas 
        write(stdout,* ) 'Tcrit                                    ',Tcrit 
        write(stdout,* ) 'Pcrit                                    ',Pcrit 
        write(stdout,* ) 'Vcrit                                    ',Vcrit 
        write(stdout,* ) '-------------------------'
        write(stdout,* ) 
        write(stdout,* ) 'Reduced parameters:'
        write(stdout,* ) '-------------------------'
        write(stdout,* ) 'Reduced Temperature      [-]             ',Tref
        write(stdout,* ) 'Reduced Pressure         [-]             ',Pref
        write(stdout,* ) 'Reduced Density          [-]             ',Rhoref
        write(stdout,* ) 'Specific gas constant    [-]             ',Rref
        write(stdout,* ) 'Heat capacity cp         [-]             ',Cpref
        write(stdout,* ) 'Speed of sound           [-]             ',SOSref
      end select
        write(stdout,* ) '-------------------------'
        write(stdout,* )  

      ! Viscosity and conductivity
      write(stdout,* )   
      write(stdout,* ) 'Viscosity and conductivity parameters:' 
      write(stdout,* ) 'Viscosity and conductivity law:        ', USE_VISC
      select case (USE_VISC)
        case("Sutherland")  ! polytropic ideal gas
          write(stdout,* ) 'Ref. temperature         [K]             ',Stref 
          write(stdout,* ) 'Ref. viscosity *10^5     [kg/(m*s)]      ',Muref * 1.0e5_mytype
          write(stdout,* ) 'Sutherland constant visc.[K]             ',Smuref
          write(stdout,* )
          write(stdout,* ) 'Ref. conductivity        [J/(m*s*K)]     ',Kref * 1.0e5_mytype
          write(stdout,* ) 'Sutherland constant cond.[K]             ',Skref
          write(stdout,* )
        case("JST")  ! cubic gas  
          write(stdout,* )  
          write(stdout,* ) 'Ref. viscosity              [-]      ',Muref
          write(stdout,* ) 'Ref. conductivity           [-]      ',Kref 
        case("Chung")  ! cubic gas  
          write(stdout,* )  
          write(stdout,* ) 'Ref. viscosity              [Pas]      ',Muref
          write(stdout,* ) 'Ref. conductivity           [W/(mK)]   ',Kref 
      end select
      write(stdout,* ) '-------------------------'
    endif
  end subroutine


end module mod_param







