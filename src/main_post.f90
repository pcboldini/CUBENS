
!--------------------------------------------------------------------------------
!
!     POST PROCESSING CODE
!
!--------------------------------------------------------------------------------


program DNS_POST

  use decomp_2d
  use decomp_2d_fft
  use decomp_2d_io
  use decomp_2d_constants
  use factor
  use mpi

  use mod_param
  use mod_halo
  use mod_grid
  use mod_init
  ! use mod_init
  use mod_eos
  use mod_eos_var
  use mod_eos_visc
  use mod_finitediff
  use mod_rhs
  use mod_auxl
  use mod_math
  use mod_postpro
  use mod_boundary
  use iso_fortran_env

  implicit none

  integer :: fh
  integer(kind=MPI_OFFSET_KIND) :: filesize, disp
  integer :: ierr,i,j,k,ii_index
  integer :: istep,nfiles,count,size_spany
  TYPE (DECOMP_INFO) :: part1,partfx,partfz,partfy,partftime !!!
  TYPE (DECOMP_INFO) :: partinterp
  character*7 :: cha
  character*1 :: cha2
  integer, dimension(3) :: fft_start, fft_end, fft_size
  character(len=30) :: date
  real(8) :: wt_start

  ! SCRINS Version number
   real(mytype), parameter                    :: version = 1.1                      ! SCRINS version

  ! underlying FFT library only needs to be initialised once
  !logical, save :: fft_initialised = .false.

  real(mytype) :: factAvg, time,theta, dTx

  real(mytype), allocatable, dimension(:,:,:) :: rho,u,v,w,ien,pre,tem,mu,ka,Cp
  real(mytype), allocatable, dimension(:,:,:) :: qvort,omx,omy,omz,vortx,vorty,vortz
  real(mytype), allocatable, dimension(:,:,:) :: dilla2,sxx,sxy,sxz,syy,syz,szz
  real(mytype), allocatable, dimension(:,:,:) :: tmp_x_arr,tmp_y_arr, tmp_z_arr
  real(mytype), allocatable, dimension(:,:,:) :: tmpPlane, inputfft

  real(mytype), allocatable, dimension(:,:) ::  arho, aT, aP, amu, aka, arT, aCp, adwdx
  real(mytype), allocatable, dimension(:,:) ::  w_BL, p_BL, rw_BL
  real(mytype), allocatable, dimension(:,:,:) :: au, aru, auu, aruu, aTauij, aqj
  real(mytype), allocatable, dimension(:,:,:,:) :: arho_time, au_time, av_time, aw_time, aT_time, aP_time
  real(mytype), allocatable, dimension(:,:,:,:) :: aCp_time, amu_time, aka_time, arw_time
  real(mytype), allocatable, dimension(:,:,:,:) :: arho_fluc, au_fluc, av_fluc, aw_fluc, aT_fluc, aP_fluc
  real(mytype), allocatable, dimension(:,:,:,:) :: aCp_fluc, amu_fluc, aka_fluc, arw_fluc
  real(mytype), allocatable, dimension(:,:,:,:) :: aw_fluc_FFT, aP_fluc_FFT, arw_fluc_FFT
  real(mytype), allocatable, dimension(:,:,:,:) :: Tauxz,qx

  complex(mytype), allocatable, dimension(:,:,:) :: spec_rho, spec_w, spec_u, spec_rw, spec_p  !!!!!
  real(mytype), allocatable, dimension(:) :: tauw, qw, rhow, muw, kaw, Tw, theta_vector, Re_tauw
  real(mytype), allocatable, dimension(:) :: Tw_global, rhow_global, muw_global, kaw_global, tauw_global, qw_global, &
                                             Re_tauw_global, theta_global

  real(mytype), allocatable, dimension(:,:) :: u_tau, Re_tau
  real(mytype), allocatable, dimension(:,:) :: u_tau_sl, Re_tau_sl
  real(mytype), allocatable, dimension(:,:) :: tTauxz_w, tqx_w

  real(mytype), allocatable, dimension(:,:,:) ::  Prod_wfluc_wfluc, Prod_wfluc_ufluc, Prod_ufluc_ufluc, Prod_vfluc_vfluc
  real(mytype), allocatable, dimension(:,:,:) ::  Prod_rfluc_rfluc, Prod_mufluc_mufluc, Prod_kafluc_kafluc, Prod_Cpfluc_Cpfluc
  real(mytype), allocatable, dimension(:,:) ::  aProd_mufluc_mufluc, aProd_rfluc_rfluc, aProd_kafluc_kafluc, aProd_Cpfluc_Cpfluc
  real(mytype), allocatable, dimension(:,:) ::  aProd_wfluc_wfluc, aProd_wfluc_ufluc, aProd_ufluc_ufluc, aProd_vfluc_vfluc
  
  complex(mytype), allocatable, dimension(:,:,:) :: spec_rho_global, spec_w_global, spec_u_global !!!!!
  complex(mytype), allocatable, dimension(:,:,:) :: tmp_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: kk = 6
  real(mytype) :: cff
  real(mytype), dimension(kk) :: in1, out1
  real(mytype), dimension(kk) :: test_in, test_out, test_out1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

interface
   subroutine assemble_globalz1D(local1D,global1D,nz)
     use decomp_2d
     integer, intent(IN) :: nz
     real(mytype), dimension(:), intent(IN) :: local1D !!!!
     real(mytype), dimension(nz_global), intent(OUT) :: global1D !!!!
   end subroutine assemble_globalz1D
end interface


  !real(mytype) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz
  !real(mytype) :: oxx, oxy, oxz, oyy, oyz, ozz

!=====================================!
!                                     !
!      INITIALIZATION                 !
!                                     !
!=====================================!


  call fdate(date) 

!-------------------------------------!
!  Read the config and init_BL files  !
!-------------------------------------!

! Reading parameters
  call read_config()
  select case (CASE)
    case("BoundaryLayer")
      call read_initBL_params()
  end select

 nfiles = (iend_pp-istart_pp)/istep_pp+1

! init MPI and decomp_2d
  call mpi_init(ierr)
  call decomp_2d_init(imax,jmax,kmax,p_row_pp,p_col_pp)
  call splitcomm(nrank,p_row_pp,p_col_pp) ! DECOMP_2D_COMM_CART_X
  call comm_init(nrank,p_row_pp,p_col_pp,DECOMP_2D_COMM_CART_X,xsize,imax,jmax,kmax)
  call decomp_2d_fft_init(PHYSICAL_IN_X)
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

  call get_decomp_info(part1)
  call get_decomp_info(partinterp)

  
!-----------------------------!
! Print a welcome message     !
!-----------------------------!

if (nrank == 0) then

  write (stdout, *)
  write (stdout, *) "o----------------------------------------------------------------------------------o"
  write (stdout, *) "|           ________  ________  ________  ___  ________   ________                 |"
  write (stdout, *) "|          |\   ____\|\   ____\|\   __  \|\  \|\   ___  \|\   ____\                |"
  write (stdout, *) "|          \ \  \___|\ \  \___|\ \  \|\  \ \  \ \  \\ \  \ \  \___|                |"
  write (stdout, *) "|           \ \_____  \ \  \    \ \   _  _\ \  \ \  \\ \  \ \_____  \              |"
  write (stdout, *) "|            \|____|\  \ \  \____\ \  \\  \\ \  \ \  \\ \  \|____|\  \             |"
  write (stdout, *) "|              ____\_\  \ \_______\ \__\\ _\\ \__\ \__\\ \__\____\_\  \            |"
  write (stdout, *) "|             |\_________\|_______|\|__|\|__|\|__|\|__| \|__|\_________\           |"
  write (stdout, *) "|             \|_________|                                  \|_________|           |"
  write (stdout, *) "|                                                                                  |"
  write (stdout, *) "|                                                                                  |"
  write (stdout, *) "|                        SuperCRItical Navier-Stokes solver                        |"
  write (stdout, *) "|                                                                                  |"
  write (stdout, '(A,F4.2,A,A,A)') " | Version ",version,": ",date,"                                     |"                                                                         
  write (stdout, *) "|                                                                                  |" 
  write (stdout, *) "o----------------------------------------------------------------------------------o"
  write (stdout, *)
  write (stdout, '(A,A)') " compiled with ", compiler_version()
  write (stdout, *)
  write (stdout, '(A,I0,A,I0,A,I0,A)') " Using ", nproc, " MPI processes, ",p_col_pp," number of procs in z and ",p_row_pp," in y"
  write (stdout, *)
  write (stdout, *) "o----------------------------------------------------------------------------------o"
  write (stdout, *)
  write (stdout, *) "|                                     POST                                         |"
  write (stdout, *)
  write (stdout, *) "o----------------------------------------------------------------------------------o"
  write (stdout, *)

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (nrank==0) write(stdout,*) 'fft_size=', fft_size

 call decomp_info_init(nfiles/2+1,imax,kmax,partftime)
 call decomp_info_init(imax,jmax,kmax/2+1,partfz)
 call decomp_info_init(imax,jmax/2+1,kmax,partfy)
 call decomp_info_init(imax/2+1,jmax,kmax,partfx)

if (nrank==0) write(stdout,*) 'xpencil (partftime%xsz)', partfy%xsz
if (nrank==0) write(stdout,*) 'ypencil (partftime%ysz)', partfy%ysz
if (nrank==0) write(stdout,*) 'zpencil (partftime%zsz)', partfy%zsz

if (fft_flag==1)  then 
  if (nrank==0) write(stdout,*) 'index_fft_span', index_fft_span
  size_spany=size(index_fft_span)
else
  index_fft_span=1
  size_spany=1
endif

! allocate main variables 
  allocate(rho(1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(u  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(v  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(w  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(ien (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(pre(1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(tem(1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(mu (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(ka (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))

  allocate( Tw(1:partfy%ysz(3))    )
  allocate( rhow(1:partfy%ysz(3))    )
  allocate( muw(1:partfy%ysz(3))    )
  allocate( kaw(1:partfy%ysz(3))    )
  allocate( tauw(1:partfy%ysz(3))    )
  allocate( qw(1:partfy%ysz(3))    )
  allocate( Re_tauw(1:partfy%ysz(3))    )
  allocate( theta_vector(1:partfy%ysz(3)) )

  allocate( Tw_global(nz_global)   )
  allocate( rhow_global(nz_global)   )
  allocate( muw_global(nz_global)   )
  allocate( kaw_global(nz_global)   )
  allocate( tauw_global(nz_global)   )
  allocate( qw_global(nz_global)   )
  allocate( Re_tauw_global(nz_global)   )
  allocate( theta_global(nz_global) )

  allocate(Re_tau(xsize(1),xsize(3)));     Re_tau=0.0_mytype;
  allocate(Re_tau_sl(xsize(1),xsize(3)));  Re_tau_sl=0.0_mytype;
  allocate(u_tau(xsize(1),xsize(3)));      u_tau = 0.0_mytype;
  allocate(u_tau_sl(xsize(1),xsize(3)));   u_tau_sl = 0.0_mytype;

  allocate(tTauxz_w(xsize(2),xsize(3)));     tTauxz_w=0.0_mytype;
  allocate(tqx_w(xsize(2),xsize(3)));  tqx_w=0.0_mytype;

  allocate(Cp(xsize(1),xsize(2),xsize(3)))
  allocate(qvort(xsize(1),xsize(2),xsize(3))); qvort = 0.0_mytype;
  allocate(dilla2(xsize(1),xsize(2),xsize(3)))
  allocate(sxx(xsize(1),xsize(2),xsize(3)))
  allocate(sxy(xsize(1),xsize(2),xsize(3)))
  allocate(sxz(xsize(1),xsize(2),xsize(3)))
  allocate(syy(xsize(1),xsize(2),xsize(3)))
  allocate(syz(xsize(1),xsize(2),xsize(3)))
  allocate(szz(xsize(1),xsize(2),xsize(3)))

  allocate(tmp_x_arr(xsize(1),xsize(2),xsize(3)))
  allocate(tmp_y_arr(xsize(1),xsize(2),xsize(3)))
  allocate(tmp_z_arr(xsize(1),xsize(2),xsize(3)))

  allocate(omx(xsize(1),xsize(2),xsize(3)))
  allocate(omy(xsize(1),xsize(2),xsize(3)))
  allocate(omz(xsize(1),xsize(2),xsize(3)))

  allocate(tmpPlane (xsize(1),xsize(2),xsize(3)))
  allocate(tmp_complex (xsize(1),xsize(2),xsize(3)))

  allocate(    arho(  xsize(1),xsize(3)));     arho = 0.0_mytype;
  allocate(    aT(  xsize(1),xsize(3)));       aT = 0.0_mytype;
  allocate(    aP(  xsize(1),xsize(3)));       aP = 0.0_mytype;
  allocate(   amu(  xsize(1),xsize(3)));       amu = 0.0_mytype;
  allocate(   aka(  xsize(1),xsize(3)));       aka = 0.0_mytype;
  allocate(   aCp(  xsize(1),xsize(3)));       aCp = 0.0_mytype;
  allocate(   adwdx(  xsize(1),xsize(3)));     adwdx = 0.0_mytype;

  allocate(    w_BL(  xsize(1),xsize(3)));     w_BL = 0.0_mytype;
  allocate(    p_BL(  xsize(1),xsize(3)));     p_BL = 0.0_mytype;
  allocate(    rw_BL(  xsize(1),xsize(3)));    rw_BL = 0.0_mytype;

  allocate(   arT(  xsize(1),xsize(3)));       arT = 0.0_mytype;
  allocate(    au(3,xsize(1),xsize(3)));       au = 0.0_mytype; !3 components
  allocate(   aru(3,xsize(1),xsize(3)));       aru = 0.0_mytype; !3 components 
  allocate(   auu(6,xsize(1),xsize(3)));       auu = 0.0_mytype; !6 components
  allocate(  aruu(6,xsize(1),xsize(3)));       aruu = 0.0_mytype; !6 components
  allocate(aTauij(6,xsize(1),xsize(3)));       aTauij = 0.0_mytype; !6 components
  allocate(   aqj(3,xsize(1),xsize(3)));       aqj = 0.0_mytype; !3 components

  allocate(    Tauxz( nfiles,xsize(1),xsize(2),xsize(3)));        Tauxz = 0.0_mytype;
  allocate(    qx( nfiles,xsize(1),xsize(2),xsize(3)));        qx = 0.0_mytype;

  if ((fft_flag==1) .OR. (rms_flag==1))  then 

  allocate(    arho_time( nfiles,xsize(1),xsize(2),xsize(3)));      arho_time = 0.0_mytype;  
  allocate(    au_time( nfiles,xsize(1),xsize(2),xsize(3)));        au_time = 0.0_mytype;
  allocate(    av_time( nfiles,xsize(1),xsize(2),xsize(3)));        av_time = 0.0_mytype;
  allocate(    aw_time( nfiles,xsize(1),xsize(2),xsize(3)));        aw_time = 0.0_mytype;
  allocate(    aT_time( nfiles,xsize(1),xsize(2),xsize(3)));        aT_time = 0.0_mytype;
  allocate(    aP_time( nfiles,xsize(1),xsize(2),xsize(3)));        aP_time = 0.0_mytype;
  allocate(    amu_time( nfiles,xsize(1),xsize(2),xsize(3)));        amu_time = 0.0_mytype;
  allocate(    aka_time( nfiles,xsize(1),xsize(2),xsize(3)));        aka_time = 0.0_mytype;
  allocate(    aCp_time( nfiles,xsize(1),xsize(2),xsize(3)));        aCp_time = 0.0_mytype;

  allocate(    arw_time( nfiles,xsize(1),xsize(2),xsize(3)));        arw_time = 0.0_mytype;

  allocate(    arho_fluc(nfiles,xsize(1),xsize(2),xsize(3)));       arho_fluc = 0.0_mytype;
  allocate(    au_fluc(nfiles,xsize(1),xsize(2),xsize(3)));         au_fluc = 0.0_mytype;
  allocate(    av_fluc(nfiles,xsize(1),xsize(2),xsize(3)));         av_fluc = 0.0_mytype;
  allocate(    aw_fluc(nfiles,xsize(1),xsize(2),xsize(3)));         aw_fluc = 0.0_mytype;
  allocate(    aT_fluc(nfiles,xsize(1),xsize(2),xsize(3)));         aT_fluc = 0.0_mytype;
  allocate(    amu_fluc(nfiles,xsize(1),xsize(2),xsize(3)));         amu_fluc = 0.0_mytype;
  allocate(    aka_fluc(nfiles,xsize(1),xsize(2),xsize(3)));         aka_fluc = 0.0_mytype;
  allocate(    aP_fluc(nfiles,xsize(1),xsize(2),xsize(3)));         aP_fluc = 0.0_mytype;
  allocate(    aCp_fluc(nfiles,xsize(1),xsize(2),xsize(3)));         aCp_fluc = 0.0_mytype;

  allocate(    arw_fluc( nfiles,xsize(1),xsize(2),xsize(3)));        arw_time = 0.0_mytype;

  allocate(    aw_fluc_FFT(nfiles,xsize(1),xsize(2),xsize(3)));         aw_fluc_FFT = 0.0_mytype;
  allocate(    arw_fluc_FFT(nfiles,xsize(1),xsize(2),xsize(3)));        arw_fluc_FFT = 0.0_mytype;
  allocate(    aP_fluc_FFT(nfiles,xsize(1),xsize(2),xsize(3)));         aP_fluc_FFT = 0.0_mytype;

  allocate(   Prod_wfluc_wfluc(  xsize(1),xsize(2),xsize(3)));       Prod_wfluc_wfluc = 0.0_mytype;
  allocate(   Prod_wfluc_ufluc(  xsize(1),xsize(2),xsize(3)));       Prod_wfluc_ufluc = 0.0_mytype;
  allocate(   Prod_ufluc_ufluc(  xsize(1),xsize(2),xsize(3)));       Prod_ufluc_ufluc = 0.0_mytype;
  allocate(   Prod_vfluc_vfluc(  xsize(1),xsize(2),xsize(3)));       Prod_vfluc_vfluc = 0.0_mytype;
  allocate(   Prod_rfluc_rfluc(  xsize(1),xsize(2),xsize(3)));       Prod_rfluc_rfluc = 0.0_mytype;
  allocate(   Prod_mufluc_mufluc(  xsize(1),xsize(2),xsize(3)));     Prod_mufluc_mufluc = 0.0_mytype;
  allocate(   Prod_kafluc_kafluc(  xsize(1),xsize(2),xsize(3)));     Prod_kafluc_kafluc = 0.0_mytype;
  allocate(   Prod_Cpfluc_Cpfluc(  xsize(1),xsize(2),xsize(3)));     Prod_Cpfluc_Cpfluc = 0.0_mytype;

  allocate(   aProd_wfluc_wfluc(  xsize(1),xsize(3)));       aProd_wfluc_wfluc = 0.0_mytype;
  allocate(   aProd_wfluc_ufluc(  xsize(1),xsize(3)));       aProd_wfluc_ufluc = 0.0_mytype;
  allocate(   aProd_ufluc_ufluc(  xsize(1),xsize(3)));       aProd_ufluc_ufluc = 0.0_mytype;
  allocate(   aProd_vfluc_vfluc(  xsize(1),xsize(3)));       aProd_vfluc_vfluc = 0.0_mytype;
  allocate(   aProd_rfluc_rfluc(  xsize(1),xsize(3)));       aProd_rfluc_rfluc = 0.0_mytype;
  allocate(   aProd_mufluc_mufluc(  xsize(1),xsize(3)));     aProd_mufluc_mufluc = 0.0_mytype;
  allocate(   aProd_kafluc_kafluc(  xsize(1),xsize(3)));     aProd_kafluc_kafluc = 0.0_mytype;
  allocate(   aProd_Cpfluc_Cpfluc(  xsize(1),xsize(3)));     aProd_Cpfluc_Cpfluc = 0.0_mytype;

  endif

  if (fft_flag==1)  then 

  allocate( inputfft(1:partfy%xsz(1),1:partfy%xsz(2),1:partfy%xsz(3))   )
  allocate( spec_rho(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )
  allocate( spec_w(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )
  !allocate( spec_u(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )
  allocate( spec_rw(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )
  allocate( spec_p(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )

  endif


  ! initialize EOS model
  if (nrank==0) then
    write(stdout,*) "Initializing EoS"
    write(stdout,*)
  endif
  call init_EOSModel()

  ! initialize visc model
  if (nrank==0) then
    write(stdout,*) "Initializing viscosity and conductivity"
    write(stdout,*)
  endif
  call init_VISCModel()

! initialize EOS and VISC
  if (nrank==0) then
    write(stdout,*) 'Initializing EOS and model for transport coefficients'
    write(stdout,*)
  endif
  call init_PARAM_EOS()

  if (nrank==0) write(stdout,*) 'initializing grid'
  call init_grid()

   if (nrank==0) write(stdout,*) 'initializing finite difference coeffs'
  call init_derivCoeffs()

  if (nrank==0) write(stdout,*) 'initializing rhs'
  call init_rhs()

  if (nrank==0) write(stdout,*) 'p_row_pp', p_row_pp
  if (nrank==0) write(stdout,*) 'p_col_pp', p_col_pp

  if (nrank==0) write(stdout,*) 'Initializing sponge'
  call initBlasius1D(part1,x,z,rho,u,v,w,ien,pre,tem,mu,ka)

  w_BL(1:xsize(1),1:xsize(3)) = w(1:xsize(1),1,1:xsize(3))
  p_BL(1:xsize(1),1:xsize(3)) = pre(1:xsize(1),1,1:xsize(3))
  rw_BL(1:xsize(1),1:xsize(3)) = rho(1:xsize(1),1,1:xsize(3))*w(1:xsize(1),1,1:xsize(3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------
! SPANWISE AVERAGING 
!----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------------
!
! all steps
!
  count=0
  wt_start = MPI_WTIME()

   nfiles = (iend_pp-istart_pp)/istep_pp+1
   if (nrank==0) write(*,*) 'number of files to average: ', nfiles-1

  factAvg = 1.0*xsize(2)*(nfiles-1)

  !Main do loop (First)

   do istep = istart_pp, (iend_pp-istep_pp), istep_pp

     count=count+1

     write(cha,'(I0.7)') istep
     if (nrank==0) write(*,'(/,A,A)') 'reading restart: ', cha

     call loadRestart(istep,time,rho,u,v,w,ien,nHalo,partinterp)

     if (nrank==0) write(stdout,*) 'updating ien,pre,tem,mu,ka'
     call calcState_re(rho,ien,pre,tem,mu,ka, 1,xsize(1),1,xsize(2),1,xsize(3))

     if (nrank==0) write(stdout,*) 'set BC'

     call haloUpdateMult_CPU((/.false.,.true.,.true./),xsize,rho,u,v,w,ien,pre,tem,mu,ka) 
     call setMultipleHaloFor1stDerBC_X(rho,u,v,w,ien,pre,tem,mu,ka)
     
     if ((perBC(3) .eqv. .false.) .and. (neigh%inlet )) call setMultipleHaloFor1stDerBC_Z_Inlet(rho,u,v,w,ien,pre,tem,mu,ka)
     if ((perBC(3) .eqv. .false.) .and. (neigh%outlet)) call setMultipleHaloFor1stDerBC_Z_Outlet(rho,u,v,w,ien,pre,tem,mu,ka)

     if (nrank==0) write(stdout,*) 'calc q-criterion'
     call calcQ(qvort,u,v,w,part1) ! q-criterion
     call decomp_2d_write_one(1,qvort, 'postproc/results/vort/qvort.'//cha//'.bin')

     if (nrank==0) write(*,*) 'calc Cp'
     call calcCp(Cp,rho,ien)

     !if (nrank==0) write(*,*) 'calc vorticity' 
     !call calcVort(vortx,vorty,vortz,u,v,w) ! vorticity
     !call decomp_2d_write_one(1,vortz,'postproc/results/vortz.'//cha//'.bin')
     if (nrank==0) write(*,*) 'calc gradients' 
     call calcGrad(part1,istep,2,yi_plane,rho,'ypl.','gradR.') ! normalised gradient

     call calcStrain(dilla2,sxx,sxy,sxz,syy,syz,szz,u,v,w) ! stress-tensor
     call calcTemp(tmp_x_arr,tmp_y_arr,tmp_z_arr,tem) ! temperature gradient
     
     do k=1,xsize(3)
      do i=1,xsize(1)
              arho(i,k) =    arho(i,k) + sum(rho(i,1:xsize(2),k))/factAvg
              aT(i,k)   =    aT(i,k) + sum(tem(i,1:xsize(2),k))/factAvg
              aP(i,k)   =    aP(i,k) + sum(pre(i,1:xsize(2),k))/factAvg
              amu(i,k)  =    amu(i,k) + sum(mu(i,1:xsize(2),k))/factAvg
              aka(i,k)  =    aka(i,k) + sum(ka(i,1:xsize(2),k))/factAvg
              adwdx(i,k)=    adwdx(i,k) + sum(sxz(i,1:xsize(2),k))/factAvg

              au(1,i,k) =   au(1,i,k) + sum(u(i,1:xsize(2),k))/factAvg
              au(2,i,k) =   au(2,i,k) + sum(v(i,1:xsize(2),k))/factAvg
              au(3,i,k) =   au(3,i,k) + sum(w(i,1:xsize(2),k))/factAvg   

              aru(1,i,k) =  aru(1,i,k) + sum(rho(i,1:xsize(2),k)*u(i,1:xsize(2),k))/factAvg
              aru(2,i,k) =  aru(2,i,k) + sum(rho(i,1:xsize(2),k)*v(i,1:xsize(2),k))/factAvg
              aru(3,i,k) =  aru(3,i,k) + sum(rho(i,1:xsize(2),k)*w(i,1:xsize(2),k))/factAvg 

              arT(i,k)   =  arT(i,k)   + sum(rho(i,1:xsize(2),k)*tem(i,1:xsize(2),k))/factAvg        

              auu(1,i,k) =  auu(1,i,k) + sum(u(i,1:xsize(2),k)*u(i,1:xsize(2),k))/factAvg
              auu(2,i,k) =  auu(2,i,k) + sum(u(i,1:xsize(2),k)*v(i,1:xsize(2),k))/factAvg
              auu(3,i,k) =  auu(3,i,k) + sum(u(i,1:xsize(2),k)*w(i,1:xsize(2),k))/factAvg
              auu(4,i,k) =  auu(4,i,k) + sum(v(i,1:xsize(2),k)*v(i,1:xsize(2),k))/factAvg
              auu(5,i,k) =  auu(5,i,k) + sum(v(i,1:xsize(2),k)*w(i,1:xsize(2),k))/factAvg
              auu(6,i,k) =  auu(6,i,k) + sum(w(i,1:xsize(2),k)*w(i,1:xsize(2),k))/factAvg

              aruu(1,i,k) = aruu(1,i,k) + sum(rho(i,1:xsize(2),k)*u(i,1:xsize(2),k)*u(i,1:xsize(2),k))/factAvg
              aruu(2,i,k) = aruu(2,i,k) + sum(rho(i,1:xsize(2),k)*u(i,1:xsize(2),k)*v(i,1:xsize(2),k))/factAvg
              aruu(3,i,k) = aruu(3,i,k) + sum(rho(i,1:xsize(2),k)*u(i,1:xsize(2),k)*w(i,1:xsize(2),k))/factAvg
              aruu(4,i,k) = aruu(4,i,k) + sum(rho(i,1:xsize(2),k)*v(i,1:xsize(2),k)*v(i,1:xsize(2),k))/factAvg
              aruu(5,i,k) = aruu(5,i,k) + sum(rho(i,1:xsize(2),k)*v(i,1:xsize(2),k)*w(i,1:xsize(2),k))/factAvg
              aruu(6,i,k) = aruu(6,i,k) + sum(rho(i,1:xsize(2),k)*w(i,1:xsize(2),k)*w(i,1:xsize(2),k))/factAvg

              aTauij(1,i,k) = aTauij(1,i,k) + sum(mu(i,1:xsize(2),k)*sxx(i,1:xsize(2),k))/factAvg ! xx
              aTauij(2,i,k) = aTauij(2,i,k) + sum(mu(i,1:xsize(2),k)*sxy(i,1:xsize(2),k))/factAvg ! xy
              aTauij(3,i,k) = aTauij(3,i,k) + sum(mu(i,1:xsize(2),k)*sxz(i,1:xsize(2),k))/factAvg ! xz
              aTauij(4,i,k) = aTauij(4,i,k) + sum(mu(i,1:xsize(2),k)*syy(i,1:xsize(2),k))/factAvg ! yy
              aTauij(5,i,k) = aTauij(5,i,k) + sum(mu(i,1:xsize(2),k)*syz(i,1:xsize(2),k))/factAvg ! yz
              aTauij(6,i,k) = aTauij(6,i,k) + sum(mu(i,1:xsize(2),k)*szz(i,1:xsize(2),k))/factAvg ! zz

              aqj(1,i,k) = aqj(1,i,k) + sum(ka(i,1:xsize(2),k)*tmp_x_arr(i,1:xsize(2),k))/factAvg ! x
              aqj(2,i,k) = aqj(2,i,k) + sum(ka(i,1:xsize(2),k)*tmp_y_arr(i,1:xsize(2),k))/factAvg ! y
              aqj(3,i,k) = aqj(3,i,k) + sum(ka(i,1:xsize(2),k)*tmp_z_arr(i,1:xsize(2),k))/factAvg ! z

              aCp(i,k) = aCp(i,k) + sum(Cp(i,1:xsize(2),k))/factAvg

              if ((fft_flag==1) .OR. (rms_flag==1))  then 
                arho_time(count,i,1:xsize(2),k)=rho(i,1:xsize(2),k)
                au_time(count,i,1:xsize(2),k)=u(i,1:xsize(2),k)
                av_time(count,i,1:xsize(2),k)=v(i,1:xsize(2),k)
                aw_time(count,i,1:xsize(2),k)=w(i,1:xsize(2),k)
                aT_time(count,i,1:xsize(2),k)=tem(i,1:xsize(2),k)
                aP_time(count,i,1:xsize(2),k)=pre(i,1:xsize(2),k)
                amu_time(count,i,1:xsize(2),k)=mu(i,1:xsize(2),k)
                aka_time(count,i,1:xsize(2),k)=ka(i,1:xsize(2),k)
                aCp_time(count,i,1:xsize(2),k)=Cp(i,1:xsize(2),k)

                arw_time(count,i,1:xsize(2),k)=rho(i,1:xsize(2),k)*w(i,1:xsize(2),k)
              endif

              Tauxz(count,i,1:xsize(2),k)=mu(i,1:xsize(2),k)*sxz(i,1:xsize(2),k)
              qx(count,i,1:xsize(2),k)=ka(i,1:xsize(2),k)*tmp_x_arr(i,1:xsize(2),k)
       enddo
     enddo

   tmpPlane(1,:,:) = Tauxz(count,1,:,:)
   call decomp_2d_write_plane(1,tmpPlane,1,1,'.','postproc/planes/Tauxz.'//cha//'.bin','dummy')
   tmpPlane(1,:,:) = qx(count,1,:,:)
   call decomp_2d_write_plane(1,tmpPlane,1,1,'.','postproc/planes/qx.'//cha//'.bin','dummy')

   enddo




! Write out 2D-planes

  tmpPlane(:,1,:) = arho;       call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_r.bin','dummy')
  tmpPlane(:,1,:) = aT;         call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_T.bin','dummy')
  tmpPlane(:,1,:) = aP;         call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_p.bin','dummy')
  tmpPlane(:,1,:) = amu;        call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_mu.bin','dummy')
  tmpPlane(:,1,:) = aka;        call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_ka.bin','dummy')
  tmpPlane(:,1,:) = aCp;        call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_Cp.bin','dummy')

  tmpPlane(:,1,:) = adwdx;      call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_dwdx.bin','dummy')

  tmpPlane(:,1,:) = au(1,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_u.bin','dummy')
  tmpPlane(:,1,:) = au(2,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_v.bin','dummy')
  tmpPlane(:,1,:) = au(3,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_w.bin','dummy')

  tmpPlane(:,1,:) = auu(1,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_uu.bin','dummy')
  tmpPlane(:,1,:) = auu(2,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_uv.bin','dummy')
  tmpPlane(:,1,:) = auu(3,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_uw.bin','dummy')
  tmpPlane(:,1,:) = auu(4,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_vv.bin','dummy')
  tmpPlane(:,1,:) = auu(5,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_vw.bin','dummy')
  tmpPlane(:,1,:) = auu(6,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_ww.bin','dummy')

  tmpPlane(:,1,:) = aruu(1,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_ruu.bin','dummy')
  tmpPlane(:,1,:) = aruu(2,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_ruv.bin','dummy')
  tmpPlane(:,1,:) = aruu(3,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_ruw.bin','dummy')
  tmpPlane(:,1,:) = aruu(4,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_rvv.bin','dummy')
  tmpPlane(:,1,:) = aruu(5,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_rvw.bin','dummy')
  tmpPlane(:,1,:) = aruu(6,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_rww.bin','dummy')

  tmpPlane(:,1,:) = aTauij(1,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_T11.bin','dummy')
  tmpPlane(:,1,:) = aTauij(2,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_T12.bin','dummy')
  tmpPlane(:,1,:) = aTauij(3,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_T13.bin','dummy')
  tmpPlane(:,1,:) = aTauij(4,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_T22.bin','dummy')
  tmpPlane(:,1,:) = aTauij(5,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_T23.bin','dummy')
  tmpPlane(:,1,:) = aTauij(6,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_T33.bin','dummy')

  tmpPlane(:,1,:) = aqj(1,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_q1.bin','dummy')
  tmpPlane(:,1,:) = aqj(2,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_q2.bin','dummy')
  tmpPlane(:,1,:) = aqj(3,:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Yavg/Yave_q3.bin','dummy')




!!!! 2D
! We need the quantitites at the wall to compute certain quantities. 
! Lets create variables that determine values at the wall.

tauw = aTauij(3,1,:) ! mu*dudy
qw = aqj(1,1,:) ! absolute value
rhow = arho(1,:)
Tw = aT(1,:)
muw = amu(1,:)*Re
kaw = aka(1,:)*Re*Pra*Ec

do k=1,xsize(3) ! 
  do i=1,xsize(1)
      u_tau(i,k) = sqrt(aTauij(3,1,k)/arho(1,k)) ! viscous so only wall units
      u_tau_sl(i,k) = sqrt(aTauij(3,1,k)/arho(i,k)) ! sl hold for semi-local viscous length scale

      Re_tau(i,k) = u_tau(i,k) * arho(1,k)/amu(1,k)  !! to be multiplied by x(2) to obtain y+ 
      Re_tau_sl(i,k) = u_tau_sl(i,k) * arho(i,k)/amu(i,k)  !! to be multiplied by x(2) to obtain y*
  enddo
enddo

Re_tauw=Re_tau(1,:)

! Momentum thickness

do k=1,xsize(3)
  call trapzint(xsize(1),x,arho(:,k)*au(3,:,k)*( 1.0_mytype-au(3,:,k) ),theta)
  theta_vector(k)=theta
enddo

call assemble_globalz1D(Tw,Tw_global,partfy%ysz(3))
call assemble_globalz1D(rhow,rhow_global,partfy%ysz(3))
call assemble_globalz1D(muw,muw_global,partfy%ysz(3))
call assemble_globalz1D(kaw,kaw_global,partfy%ysz(3))
call assemble_globalz1D(tauw,tauw_global,partfy%ysz(3))
call assemble_globalz1D(qw,qw_global,partfy%ysz(3))
call assemble_globalz1D(Re_tauw,Re_tauw_global,partfy%ysz(3))
call assemble_globalz1D(theta_vector,theta_global,partfy%ysz(3))

if (nrank==0) then
     open(18,file = 'postproc/wall_prop.txt')
     write(18,'(26A71)') 'k - Tw - rhow - muw - kaw - Cf - qw - Re_tauw - theta_global'
     do k=1,nz_global
        write(18,*) k, Tw_global(k), rhow_global(k), muw_global(k), kaw_global(k), tauw_global(k), qw_global(k), &
                    Re_tauw_global(k), theta_global(k)     
     enddo
     close(18)
endif


!!!! 3D
! time average

do k=1,xsize(3)
    do j=1,xsize(2)

      tTauxz_w(j,k) = sum(Tauxz(1:count,1,j,k))/(nfiles-1)
      tqx_w(j,k) = sum(qx(1:count,1,j,k))/(nfiles-1)

  enddo
enddo

 tmpPlane(:,1,:) = Re_tau;  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Re_tau.bin','dummy')
 tmpPlane(:,1,:) = u_tau;  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/u_tau.bin','dummy')

 tmpPlane(:,1,:) = Re_tau_sl;  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/Re_tau_sl.bin','dummy')
 tmpPlane(:,1,:) = u_tau_sl;  call decomp_2d_write_plane(1,tmpPlane,2,1,'.','postproc/results/u_tau_sl.bin','dummy')

 tmpPlane(1,:,:) = tTauxz_w;  call decomp_2d_write_plane(1,tmpPlane,1,1,'.','postproc/results/tTauxz_w.bin','dummy')
 tmpPlane(1,:,:) = tqx_w;  call decomp_2d_write_plane(1,tmpPlane,1,1,'.','postproc/results/tqx_w.bin','dummy')


 if (nrank==0) write(*,*)
 if (nrank==0) write(*,*) 'Wall done!'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------
! PERTURBATIONS AND RMS
!----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (rms_flag==1)  then 

do ii_index=1,count
  do k=1,xsize(3)
    do i=1,xsize(1)
      arho_fluc(ii_index,i,1:xsize(2),k)=arho_time(ii_index,i,1:xsize(2),k)-arho(i,k)
      aka_fluc(ii_index,i,1:xsize(2),k)=aka_time(ii_index,i,1:xsize(2),k)-aka(i,k)
      amu_fluc(ii_index,i,1:xsize(2),k)=amu_time(ii_index,i,1:xsize(2),k)-amu(i,k)
      au_fluc(ii_index,i,1:xsize(2),k)=au_time(ii_index,i,1:xsize(2),k)-au(1,i,k)
      av_fluc(ii_index,i,1:xsize(2),k)=av_time(ii_index,i,1:xsize(2),k)-au(2,i,k)
      aw_fluc(ii_index,i,1:xsize(2),k)=aw_time(ii_index,i,1:xsize(2),k)-au(3,i,k)
      aP_fluc(ii_index,i,1:xsize(2),k)=aP_time(ii_index,i,1:xsize(2),k)-aP(i,k)  
      aCp_fluc(ii_index,i,1:xsize(2),k)=aCp_time(ii_index,i,1:xsize(2),k)-aCp(i,k)   

      arw_fluc(ii_index,i,1:xsize(2),k)=arw_time(ii_index,i,1:xsize(2),k)-aru(3,i,k)
                
    enddo
  enddo    
enddo 

do k =1,xsize(3)
  do i =1,xsize(1)
    do j = 1, xsize(2)
    ! Taking all the Time, FactAvg will average on  both time and Span
      do ii_index=1,count

        Prod_wfluc_wfluc(i,j,k) = Prod_wfluc_wfluc(i,j,k) +  aw_fluc(ii_index,i,j,k)*aw_fluc(ii_index,i,j,k) ! front files j = 1
        Prod_wfluc_ufluc(i,j,k) = Prod_wfluc_ufluc(i,j,k) +  aw_fluc(ii_index,i,j,k)*au_fluc(ii_index,i,j,k)
        Prod_ufluc_ufluc(i,j,k) = Prod_ufluc_ufluc(i,j,k) +  au_fluc(ii_index,i,j,k)*au_fluc(ii_index,i,j,k)
        Prod_vfluc_vfluc(i,j,k) = Prod_vfluc_vfluc(i,j,k) +  av_fluc(ii_index,i,j,k)*av_fluc(ii_index,i,j,k)
        Prod_rfluc_rfluc(i,j,k) = Prod_rfluc_rfluc(i,j,k) +  arho_fluc(ii_index,i,j,k)*arho_fluc(ii_index,i,j,k)
        Prod_mufluc_mufluc(i,j,k) = Prod_mufluc_mufluc(i,j,k) +  amu_fluc(ii_index,i,j,k)*amu_fluc(ii_index,i,j,k)
        Prod_kafluc_kafluc(i,j,k) = Prod_kafluc_kafluc(i,j,k) +  aka_fluc(ii_index,i,j,k)*aka_fluc(ii_index,i,j,k)
        Prod_Cpfluc_Cpfluc(i,j,k) = Prod_Cpfluc_Cpfluc(i,j,k) +  aCp_fluc(ii_index,i,j,k)*aCp_fluc(ii_index,i,j,k)

      enddo
    enddo
    
    ! Span and Time average
     aProd_wfluc_wfluc(i,k) = aProd_wfluc_wfluc(i,k) + sum(Prod_wfluc_wfluc(i,1:xsize(2),k))/factAvg
     aProd_wfluc_ufluc(i,k) = aProd_wfluc_ufluc(i,k) + sum(Prod_wfluc_ufluc(i,1:xsize(2),k))/factAvg
     aProd_ufluc_ufluc(i,k) = aProd_ufluc_ufluc(i,k) + sum(Prod_ufluc_ufluc(i,1:xsize(2),k))/factAvg
     aProd_vfluc_vfluc(i,k) = aProd_vfluc_vfluc(i,k) + sum(Prod_vfluc_vfluc(i,1:xsize(2),k))/factAvg
     aProd_rfluc_rfluc(i,k) = aProd_rfluc_rfluc(i,k) + sum(Prod_rfluc_rfluc(i,1:xsize(2),k))/factAvg
     aProd_mufluc_mufluc(i,k) = aProd_mufluc_mufluc(i,k)+ sum(Prod_mufluc_mufluc(i,1:xsize(2),k))/factAvg
     aProd_kafluc_kafluc(i,k) = aProd_kafluc_kafluc(i,k)+ sum(Prod_kafluc_kafluc(i,1:xsize(2),k))/factAvg
     aProd_Cpfluc_Cpfluc(i,k) = aProd_Cpfluc_Cpfluc(i,k) + sum(Prod_Cpfluc_Cpfluc(i,1:xsize(2),k))/factAvg

    ! Wall units
    ! sa_wfluc_wfluc_p(i,k) = sa_wfluc_wfluc(i,k)*u_tau(i,k)*arho(i,k)/amu(i,k) 
    ! sa_wfluc_ufluc_p(i,k) = sa_wfluc_ufluc(i,k)*u_tau(i,k)*arho(i,k)/amu(i,k) 
    ! sa_ufluc_ufluc_p(i,k) = sa_ufluc_ufluc(i,k)*u_tau(i,k)*arho(i,k)/amu(i,k) 
    ! sa_vfluc_vfluc_p(i,k) = sa_vfluc_vfluc(i,k)*u_tau(i,k)*arho(i,k)/amu(i,k) 

  enddo
enddo


tmpPlane(:,1,:) = aProd_wfluc_wfluc(:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1, &
                                                                     '.','postproc/results/rms/aProd_wfluc_wfluc.bin','dummy')
tmpPlane(:,1,:) = aProd_wfluc_ufluc(:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1, &
                                                                     '.','postproc/results/rms/aProd_wfluc_ufluc.bin','dummy')
tmpPlane(:,1,:) = aProd_ufluc_ufluc(:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1, &
                                                                     '.','postproc/results/rms/aProd_ufluc_ufluc.bin','dummy')
tmpPlane(:,1,:) = aProd_vfluc_vfluc(:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1, &
                                                                     '.','postproc/results/rms/aProd_vfluc_vfluc.bin','dummy')
tmpPlane(:,1,:) = aProd_rfluc_rfluc(:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1, &
                                                                      '.','postproc/results/rms/aProd_rfluc_rfluc.bin','dummy')
tmpPlane(:,1,:) = aProd_mufluc_mufluc(:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1, &
                                                                       '.','postproc/results/rms/aProd_mufluc_mufluc.bin','dummy')
tmpPlane(:,1,:) = aProd_kafluc_kafluc(:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1, &
                                                                       '.','postproc/results/rms/aProd_kafluc_kafluc.bin','dummy')
tmpPlane(:,1,:) = aProd_Cpfluc_Cpfluc(:,:);  call decomp_2d_write_plane(1,tmpPlane,2,1, &
                                                                       '.','postproc/results/rms/aProd_Cpfluc_Cpfluc.bin','dummy')


if (nrank==0) write(*,*) 'RMS done!'

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------
! FFT
!----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (fft_flag==1)  then 

  if (nrank == 0) then

    write(stdout,* )
    write(stdout,* ) 'FFT'
    write(stdout,* ) '-------------------------'

  endif 


  do ii_index=1,count
    do k=1,xsize(3)
      do i=1,xsize(1)
            !arho_fluc(ii_index,i,1:xsize(2),k)=arho_time(ii_index,i,1:xsize(2),k)-arho(i,k)
            !au_fluc(ii_index,i,1:xsize(2),k)=au_time(ii_index,i,1:xsize(2),k)-au(1,i,k)
            !av_fluc(ii,i,1:xsize(2),k)=av_time(ii,i,1:xsize(2),k)-au(2,i,k)
            aw_fluc_FFT(ii_index,i,1:xsize(2),k)=aw_time(ii_index,i,1:xsize(2),k)-w_BL(i,k)
            !aT_fluc(ii,i,1:xsize(2),k)=aT_time(ii,i,1:xsize(2),k)-aT(i,k)
            aP_fluc_FFT(ii_index,i,1:xsize(2),k)=aP_time(ii_index,i,1:xsize(2),k)-p_BL(i,k) 

            arw_fluc_FFT(ii_index,i,1:xsize(2),k)=arw_time(ii_index,i,1:xsize(2),k)-rw_BL(i,k)
      enddo
    enddo
    write(cha,'(I0.7)') istart_pp+istep_pp*(ii_index - 1)
    if (nrank==0) write(*,*) cha
    !call spectray(2,spec_rho,jmax,part1,partfy,arho_fluc(ii_index,1:xsize(1),1:xsize(2),1:xsize(3)))
    call spectray(2,spec_w,jmax,part1,partfy,aw_fluc_FFT(ii_index,1:xsize(1),1:xsize(2),1:xsize(3))) 
    !call spectray(2,spec_u,jmax,part1,partfy,au_fluc(ii_index,1:xsize(1),1:xsize(2),1:xsize(3)))
    call spectray(2,spec_rw,jmax,part1,partfy,arw_fluc_FFT(ii_index,1:xsize(1),1:xsize(2),1:xsize(3)))
    call spectray(2,spec_p,jmax,part1,partfy,aP_fluc_FFT(ii_index,1:xsize(1),1:xsize(2),1:xsize(3))) 


    do j=1,size_spany
      write(cha2,'(I0)') index_fft_span(j)
      ! tmp_complex = spec_rho(1:xsize(1),1:xsize(2),1:xsize(3))
      ! call decomp_2d_write_plane(1,tmp_complex,2,index_fft_span(j)+1, &
      !                           '.','postproc/results/fft/fft_r_span'//cha2//'_'//cha//'.bin', 'dummy')
      tmp_complex = spec_w(1:xsize(1),1:xsize(2),1:xsize(3))
      call decomp_2d_write_plane(1,tmp_complex,2,index_fft_span(j)+1, &
                                '.','postproc/results/fft/fft_w_span'//cha2//'_'//cha//'.bin', 'dummy')
      tmp_complex = spec_rw(1:xsize(1),1:xsize(2),1:xsize(3))
      call decomp_2d_write_plane(1,tmp_complex,2,index_fft_span(j)+1, &
                                '.','postproc/results/fft/fft_rw_span'//cha2//'_'//cha//'.bin', 'dummy')
      tmp_complex = spec_p(1:xsize(1),1:xsize(2),1:xsize(3))
      call decomp_2d_write_plane(1,tmp_complex,2,index_fft_span(j)+1, &
                                '.','postproc/results/fft/fft_p_span'//cha2//'_'//cha//'.bin', 'dummy')
      ! tmp_complex = spec_u(1:xsize(1),1:xsize(2),1:xsize(3))
      ! call decomp_2d_write_plane(1,tmp_complex,2,index_fft_span(j)+1, &
      !                           '.','postproc/results/fft/fft_u_span'//cha2//'_'//cha//'.bin', 'dummy')     
    enddo
    
  enddo 

if (nrank==0) write(*,*) 'FFT done!'

end if

if (nrank == 0) write(*,*)
if (nrank == 0) write (stdout, *) "o----------------------------------------------------------------------------------o"
if (nrank == 0) write(*,*) 'POST done!'
if (nrank == 0) print '("Total time = ",f10.3," minutes.")', (MPI_WTIME() - wt_start)/60.0

  deallocate(rho)
  deallocate(u)
  deallocate(v)
  deallocate(w)
  deallocate(ien)
  deallocate(pre)
  deallocate(tem)
  deallocate(mu)
  deallocate(ka)

  deallocate(Cp)
  deallocate(qvort)
  deallocate(dilla2)
  deallocate(sxx)
  deallocate(sxy)
  deallocate(sxz)
  deallocate(syy)
  deallocate(syz)
  deallocate(szz)

  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call mpi_finalize(ierr)

end program


subroutine assemble_globalz1D(local1D,global1D,nz)
  
  use decomp_2d
  use mpi
  
  implicit none
  
  integer, intent(IN) :: nz
  real(mytype), dimension(:), intent(IN) :: local1D !!!!
  real(mytype), dimension(nz_global), intent(OUT) :: global1D !!!!
  
  real(mytype), allocatable, dimension(:) :: rbuf_1D !!!!!
  integer, dimension(3) :: sbuf1D, rbuf1D
  
  integer :: ierror,k,m,k1,k2, count
  integer, dimension(MPI_STATUS_SIZE) :: status
  
  if (nrank==0) then
     ! master writes its own data to a global array
        k1 = 1
        k2 = nz
     do k=k1,k2
        ! 'local' is assumed shape array
        ! but it is OK as starting index for rank 0 always 1
        global1D(k)=local1D(k)
     end do
     ! then loop through all other ranks to collect data
     do m=1,nproc-1
        CALL MPI_RECV(rbuf1D,3,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
             status,ierror)
        allocate(rbuf_1D(rbuf1D(1):rbuf1D(2)))
        CALL MPI_RECV(rbuf_1D,rbuf1D(3),real_type,m, &
             m+nproc,MPI_COMM_WORLD,status,ierror) !!!!!
        do k=rbuf1D(1),rbuf1D(2)
          global1D(k)=rbuf_1D(k)
        end do
        deallocate(rbuf_1D)
     end do
  else
     ! slaves send data to master
        sbuf1D(1) = nrank*nz+1
        sbuf1D(2) = (nrank+1)*nz
        sbuf1D(3) = (nrank+1)*nz
        count = nz
     ! send partition information
     CALL MPI_SEND(sbuf1D,3,MPI_INTEGER,0,nrank,MPI_COMM_WORLD,ierror)
     ! send data array
     CALL MPI_SEND(local1D,count,real_type,0, &
          nrank+nproc,MPI_COMM_WORLD,ierror) !!!!!
  end if
  
  return
end subroutine assemble_globalz1D



