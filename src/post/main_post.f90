! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
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
  use mod_eos
  use mod_eos_var
  use mod_eos_visc
  use mod_finitediff
  use mod_rhs
  use mod_auxl
  use mod_perturbation
  use mod_math
  use mod_auxlpro
  use mod_boundary
  use iso_fortran_env
  implicit none

!===============================================================================================!
!                                       CUBENS post data
!===============================================================================================!
  integer(kind=MPI_OFFSET_KIND) :: filesize, disp
  integer :: fh, ierr,i,j,k,ii_index
  integer :: istep,nfiles,count,size_spany
  integer :: nfiles_rms,nfiles_fft,kstart,kend,kglobal,iname,lenr,izpl
  TYPE (DECOMP_INFO) :: part1,partfy 
  TYPE (DECOMP_INFO) :: partinterp
  character*7 :: cha
  character*1 :: cha2
  character(len=1024) :: cha3
  integer, dimension(3) :: fft_start, fft_end, fft_size
  character(len=30) :: date
  real(8) :: wt_start
  logical :: exist1,exist2,exist3
! CUBENS Version number
  real(mytype), parameter                    :: version = 1.0
  real(mytype) :: factAvg, time,theta, dTx
  ! 3-D properties
  real(mytype), allocatable, dimension(:,:,:) :: rho,u,v,w,ien,pre,tem,mu,ka,Cp
  real(mytype), allocatable, dimension(:,:,:) :: rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl
  real(mytype), allocatable, dimension(:,:,:) :: qvort,omx,omy,omz,vortx,vorty,vortz
  real(mytype), allocatable, dimension(:,:,:) :: dilla2,sxx,sxy,sxz,syy,syz,szz
  real(mytype), allocatable, dimension(:,:,:) :: tmp_x_arr,tmp_y_arr, tmp_z_arr
  real(mytype), allocatable, dimension(:,:,:) :: tmpPlane, inputfft
  ! spanwise averaged properties
  character(5) :: stat_name(27)
  real(mytype), allocatable, dimension(:,:) ::  arho, aT, aP, amu, aka, arT, aCp, adwdx
  real(mytype), allocatable, dimension(:,:,:) :: au, aru, auu, aruu, aTauij, aqj
  ! 3-D and time properties
  real(mytype), allocatable, dimension(:,:,:,:) :: arho_time, au_time, av_time, aw_time, aT_time, aP_time
  real(mytype), allocatable, dimension(:,:,:,:) :: aCp_time, amu_time, aka_time, arw_time
  real(mytype), allocatable, dimension(:,:,:,:) :: Tauxz_time,qx_time
  !RMS (tauw)
  real(mytype), allocatable, dimension(:)   ::  tauw_rms, tauw_rms_global
  real(mytype) :: tauw_local
  !RMS (quadrant)
  real(mytype), allocatable, dimension(:,:,:) :: aQ4
  integer, allocatable, dimension(:,:,:) :: countQ4
  real(mytype) :: u_flu, w_flu
  ! FFT properties
  real(mytype), allocatable, dimension(:,:,:) :: arho_fluc, aw_fluc, au_fluc
  real(mytype), allocatable, dimension(:,:,:,:) :: aw_fluc_FFT, aP_fluc_FFT, arw_fluc_FFT
  complex(mytype), allocatable, dimension(:,:,:) :: spec_rho, spec_w, spec_u, spec_rw, spec_p 
  complex(mytype), allocatable, dimension(:,:,:) :: tmp_complex
  ! Wall properties
  real(mytype), allocatable, dimension(:) :: tauw, qw, rhow, muw, kaw, Tw, theta_vector, Re_tauw
  real(mytype), allocatable, dimension(:) :: Tw_global, rhow_global, muw_global, kaw_global, tauw_global, qw_global, &
                                             Re_tauw_global, theta_global
  ! Turbulent properties
  real(mytype), allocatable, dimension(:,:) :: u_tau, Re_tau
  real(mytype), allocatable, dimension(:,:) :: u_tau_sl, Re_tau_sl
  real(mytype), allocatable, dimension(:,:) :: tTauxz_w, tqx_w
! Interface for global variables
interface
  subroutine assemble_globalz1D(local1D,global1D,kstart,kend)
    use decomp_2d
    integer, intent(IN) :: kstart,kend
    real(mytype), dimension(:), intent(IN) :: local1D !!!!
    real(mytype), dimension(nz_global), intent(OUT) :: global1D !!!!
  end subroutine assemble_globalz1D
end interface

!===============================================================================================!
!                                       INITIALIZATION
!===============================================================================================!
! generate today's date
  call fdate(date) 
!-------------------------------------!
!  Read the config and init_BL files  !
!-------------------------------------!
! Reading parameters
  call read_config()
#if defined(BL)
  call read_init_params()
  ! Writing parameters for user variation 
  call io_writeParams_variation()
#endif
! Calculate the number of files to average
 nfiles = (iend_pp-istart_pp)/istep_pp+1
!-----------------------------!
!         MPI library         !
!-----------------------------!
! init MPI and decomp_2d
  call mpi_init(ierr)
! p_row_pp, p_col_pp needs to be equal in config.h and when launching ./post
  call decomp_2d_init(imax,jmax,kmax,p_row_pp,p_col_pp)
  call comm_init(nrank,p_row_pp,p_col_pp,DECOMP_2D_COMM_CART_X,xsize,imax,jmax,kmax)
  call decomp_2d_fft_init(PHYSICAL_IN_X)
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
  call get_decomp_info(part1)
  call get_decomp_info(partinterp)
  call decomp_info_init(imax,jmax/2+1,kmax,partfy)
! Print a welcome message 
if (nrank == 0) then
    write (stdout, *)
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *) "|         ________  ______  _______   ________     |"
    write (stdout, *) "|        / ____/ / / / __ )/ ____/ | /  / ___/     |"
    write (stdout, *) "|       / /   / / / / /_/ / __/ /  |/  //____      |"
    write (stdout, *) "|      / /___/ /_/ / /_/ / /___/ __   /___/ /      |"
    write (stdout, *) "|     /_____/_____/_____/_____/_/  |_//____/       |"
    write (stdout, *) "|                                                  |"
    write (stdout, *) "|                                                  |"
    write (stdout, *) "|       CUBic Equation of state Navier-Stokes      |"
    write (stdout, *) "|                                                  |"
    write (stdout, '(A,F4.2,A,A,A)') " |      Version ",version,": ",date,"|"                                                                         
    write (stdout, *) "|                                                  |" 
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *)
    write (stdout, '(A,A)') " compiled with ", compiler_version()
    write (stdout, *)
    write (stdout, '(A,I0,A,I0,A,I0,A)') " Using ", nproc, " MPI processes, ", p_col_pp, &
                   " number of procs in z and ", p_row_pp, " in y " 
    write (stdout, *)
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *) "|                       POST                       |"
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *)
endif

!===============================================================================================!
!                                       ALLOCATION
!===============================================================================================!
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
  allocate(arho(  xsize(1),xsize(3)));     arho = 0.0_mytype;
  allocate(  ap(  xsize(1),xsize(3)));     ap   = 0.0_mytype;
  allocate(  aT(  xsize(1),xsize(3)));     aT   = 0.0_mytype;
  allocate( amu(  xsize(1),xsize(3)));     amu  = 0.0_mytype;
  allocate( aka(  xsize(1),xsize(3)));     aka = 0.0_mytype;
  allocate( arT(  xsize(1),xsize(3)));     arT  = 0.0_mytype;
  allocate(    au(3,xsize(1),xsize(3)));   au     = 0.0_mytype;  !3 components
  allocate(   aru(3,xsize(1),xsize(3)));   aru    = 0.0_mytype;  !3 components
  allocate(  aruu(6,xsize(1),xsize(3)));   aruu   = 0.0_mytype;  !6 components
  allocate(aTauij(6,xsize(1),xsize(3)));   aTauij = 0.0_mytype;  !6 components
  allocate(   aqj(3,xsize(1),xsize(3)));   aqj    = 0.0_mytype;  !3 components
if (avg_flag == 'restart')  then
  allocate(   aCp(  xsize(1),xsize(3)));       aCp = 0.0_mytype;
  allocate(   adwdx(  xsize(1),xsize(3)));     adwdx = 0.0_mytype;
  allocate(   auu(6,xsize(1),xsize(3)));       auu = 0.0_mytype; !6 components
  allocate(tTauxz_w(xsize(2),xsize(3)));     tTauxz_w=0.0_mytype;
  allocate(tqx_w(xsize(2),xsize(3)));        tqx_w=0.0_mytype;
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
  allocate(Tauxz_time( nfiles,xsize(1),xsize(2),xsize(3)));        Tauxz_time = 0.0_mytype;
  allocate(qx_time( nfiles,xsize(1),xsize(2),xsize(3)));        qx_time = 0.0_mytype;
endif
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

! for FFT or RMS properties
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
  allocate(    aw_fluc_FFT(nfiles,xsize(1),xsize(2),xsize(3)));         aw_fluc_FFT = 0.0_mytype;
  allocate(    arw_fluc_FFT(nfiles,xsize(1),xsize(2),xsize(3)));        arw_fluc_FFT = 0.0_mytype;
  allocate(    aP_fluc_FFT(nfiles,xsize(1),xsize(2),xsize(3)));         aP_fluc_FFT = 0.0_mytype;
  endif
! for FFT properties
  if (fft_flag==1)  then 
  allocate( inputfft(1:partfy%xsz(1),1:partfy%xsz(2),1:partfy%xsz(3))   )
  allocate( spec_rho(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )
  allocate( spec_w(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )
  allocate( spec_u(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )
  allocate( spec_rw(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )
  allocate( spec_p(1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3))   )
  endif
! Print the initial simulation parameters
  call print_init_params()
! initialize EoS 
  if (nrank==0) then
    write(stdout,*) "Initializing EoS"
  endif
  call init_EOSModel()
! initialize transport properties
  if (nrank==0) then
    write(stdout,*) "Initializing viscosity and conductivity"
  endif
  call init_VISCModel()
! initialize EoS and TP global parameters
  if (nrank==0) then
    write(stdout,*) 'Initializing global variables for EOS and transport properties'
    write(stdout,*)
  endif
  call init_PARAM_EOS()
! initialize grid
  if (nrank==0) then
    write(stdout,*) 'Initializing grid'
    write(stdout,*)
  endif
  call init_grid()
!Get partition information
  do k=1,nz_global
    if (z(1) .eq. z_global(k)) kstart = k
    if (z(xsize(3)) .eq. z_global(k)) kend = k
  enddo
  call print_init_grid()
! initialize finite-difference coefficients
  if (nrank==0) then
    write(stdout,*) 'Initializing finite-difference coefficients'
    write(stdout,*)
  endif
  call init_derivCoeffs()
  call print_init_derivCoeffs()
! allocate variables for right-hand side 
  if (nrank==0) then
    write(stdout,*) 'Initializing RHS'
    write(stdout,*)
  endif
  call init_rhs()
  call print_init_rhs() 
  ! initialize perturbation for transitional boundary layer
  if (nrank==0 .and. pert_calc==1) then
    write(stdout,*) 'Initializing perturbation'
    write(stdout,*)
  endif
  call init_pert()
! initialize initial condition for boundary layer
# if defined(BL)
  if (fft_flag==1)  then 
    if (nrank==0) then
      write(stdout,*) 'Initializing laminar BL profiles'
      write(stdout,*)
    endif
    allocate(rho_bl  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(u_bl    (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(v_bl    (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(w_bl    (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(ien_bl  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(pre_bl  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(tem_bl  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(mu_bl   (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(ka_bl   (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    call initField_BL(part1,x,z,rho,u,v,w,ien,pre,tem,mu,ka)
    ! save boundary layer profiles for FFT
    rho_bl=rho;u_bl=u;v_bl=v;w_bl=w;ien_bl=ien;pre_bl=pre;tem_bl=tem;mu_bl=mu;ka_bl=ka
  endif
#endif
! initialize restart averaging
! initialize planes averaging
if (avg_flag == 'planes') then
  if (nrank==0) then
    write(stdout, *) "o--------------------------------------------------o"
    write(stdout,*) 'Initializing planes averaging'
    write(stdout,*)
    write(stdout,'(A, I10)') 'stats_step:                           ',stats_step
    write(stdout,'(A, F10.3)') 'stats_time_rate:                      ',stats_time_rate
  endif
endif
! initialize RMS
if (rms_flag==1) then
  if (nrank==0) then
    write(stdout,*)
    write(stdout,*) 'Initializing RMS'
    write(stdout,*)
    write(stdout,'(A, I10)') 'istart_rms:                           ', istart_rms
    write(stdout,'(A, I10)') 'iend_rms:                             ', iend_rms
    write(stdout,'(A, I10)') 'index_rms_xpl:                        ', index_rms_xpl
    write(stdout,'(A, I10)') 'index_rms_zpl:                        ', index_rms_zpl
  endif
endif
! initialize FFT
  if (fft_flag==1)  then 
    if (nrank==0) then
      write(stdout,*)
      write(stdout,*) 'Initializing FFT'
      write(stdout,*)
      write (stdout, '(A,I0,A,I0,A,I0,A)') " FFT setup: ", fft_size(1), " points in x, ", fft_size(2), &
                   " points in y and ", p_row_pp, " points in z" 
    endif
    size_spany=size(index_fft_span)
  else
    index_fft_span=1
    size_spany=1
  endif
  write(stdout, *) "o--------------------------------------------------o"

!===============================================================================================!
!                                      SPANWISE AVERAGING
!===============================================================================================!
! through all steps
  count=0
  wt_start = MPI_WTIME()
  if (avg_flag == 'restart') then
    write(stdout,*) 'Average with restart files'
    write(stdout,*) 
    nfiles = (iend_pp-istart_pp)/istep_pp+1
    if (nrank==0) write(stdout,'(A, I10)') 'number of files to time average:      ', nfiles-1
    factAvg = 1.0*xsize(2)*(nfiles-1)
    ! main time loop
    do istep = istart_pp, (iend_pp-istep_pp), istep_pp
      count=count+1
      write(cha,'(I0.7)') istep
      if (nrank==0) write(stdout,'(A, A10)') 'reading restart at step:              ', cha
      call loadRestart(istep,time,rho,u,v,w,ien,nHalo,partinterp)
      if (nrank==0) write(stdout,*) 'updating ien,pre,tem,mu,ka,cp'
      call calcState_re(rho,ien,pre,tem,mu,ka, 1,xsize(1),1,xsize(2),1,xsize(3))
      call calcCp(Cp,rho,ien)
      if (nrank==0) write(stdout,*) 'setting boundary conditions'
      call setBC(part1,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
      if (nrank==0) write(stdout,*) 'calculating Q-criterion'
      ! Q-criterion
      call calcQ(qvort,u,v,w,part1) 
      call decomp_2d_write_one(1,qvort, 'postproc/results/vort/qvort.'//cha//'.bin')
      !if (nrank==0) write(*,*) 'calculating vorticity' 
      !call calcVort(vortx,vorty,vortz,strxz,u,v,w) ! vorticity
      !call decomp_2d_write_one(1,vortz,'postproc/results/vortz.'//cha//'.bin')
      if (nrank==0) write(stdout,*) 'calculate wall-normal gradients' 
      ! normalised gradient
      if (yi_plane(1).gt.0)  then
        do i=1,size(yi_plane)
          call calcGrad(part1,istep,2,yi_plane(i),rho,'ypl.','gradR.')
        enddo
      endif
      ! calculate stress-tensor
       call calcStrain(dilla2,sxx,sxy,sxz,syy,syz,szz,u,v,w) 
      ! calculate temperature gradient
       call calcTemp(tmp_x_arr,tmp_y_arr,tmp_z_arr,tem) 
      ! streamwise and wall-normal loop for spanwise averaging
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
          Tauxz_time(count,i,1:xsize(2),k)=mu(i,1:xsize(2),k)*sxz(i,1:xsize(2),k)
          qx_time(count,i,1:xsize(2),k)=ka(i,1:xsize(2),k)*tmp_x_arr(i,1:xsize(2),k)
        enddo
      enddo
      ! write planes for instantaneous wall stress and heat flux
      tmpPlane(1,:,:) = Tauxz_time(count,1,:,:)
      call decomp_2d_write_plane(1,tmpPlane,1,1,'.','postproc/planes/Tauxz_time.'//cha//'.bin','dummy')
      tmpPlane(1,:,:) = qx_time(count,1,:,:)
      call decomp_2d_write_plane(1,tmpPlane,1,1,'.','postproc/planes/qx_time.'//cha//'.bin','dummy')
    enddo

!===============================================================================================!
!                                         OUTPUT
!===============================================================================================!
! Write out 2D-planes of spanwise averaged quantities (written in /postproc/results/Yavg)
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
  ! 2-D
  do k=1,xsize(3)
    do j=1,xsize(2)
      tTauxz_w(j,k) = sum(Tauxz_time(1:count,1,j,k))/(nfiles-1)
      tqx_w(j,k) = sum(qx_time(1:count,1,j,k))/(nfiles-1)
    enddo
  enddo
  tmpPlane(1,:,:) = tTauxz_w;  call decomp_2d_write_plane(1,tmpPlane,1,1,'.','postproc/results/tTauxz_w.bin','dummy')
  tmpPlane(1,:,:) = tqx_w;  call decomp_2d_write_plane(1,tmpPlane,1,1,'.','postproc/results/tqx_w.bin','dummy')
  deallocate(tmpPlane)
  else if (avg_flag == 'planes') then
    write(stdout,*) 'Average with planes files'
    !Check number of stats
    if (nrank==0) write(stdout,'(A, I10)') 'number of stats files:                ', size(stats_step)
    if (nrank==0) write(stdout,*) 
    if (size(stats_step) .ne. size(stats_time_rate)) then
      if(nrank==0) write(stdout,*) 'Check number of stats in config.h!'
      call MPI_FILE_CLOSE(fh,ierr)
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      stop
    end if
    factAvg = 1.0_mytype / sum(stats_time_rate)
    !Set file name (number of elements are decided when the variable is defined)
    stat_name = ['r    ','u    ','v    ','w    ','p    ','t    ','mu   ', 'ka   ', &
                 'ru   ','rv   ','rw   ','rt   ',&
                 'ruu  ','ruv  ','ruw  ','rvv  ','rvw  ','rww  ', &
                 'tauxx','tauxy','tauxz','tauyy','tauyz','tauzz', &
                 'qx   ','qy   ','qz   ']
    !Allocation
    allocate(tmpPlane (xsize(1),xsize(3),size(stat_name))); tmpPlane = 0.0_mytype;
    !Main do loop
    do istep = 1,size(stats_step)
      count=count+1
      write(cha,'(I0.7)') stats_step(istep)
      if (nrank==0) write(stdout,'(A, A10)') 'reading stats:                        ', cha
      do iname = 1,size(stat_name)
        inquire(iolength=lenr) tmpPlane(1,1,iname)
        open(10,file='postproc/stats/'//trim(stat_name(iname))//'_avg.'//cha//'.bin',&
                form='unformatted',status="old",access='direct',recl=xsize(1)*lenr)
        do k=1,xsize(3)
          read(10, rec=k+kstart-1) (tmpPlane(i,k,iname), i=1,xsize(1))
        enddo
        close(10)
      enddo
      tmpPlane = tmpPlane*stats_time_rate(istep)*factAvg
      arho(:,:) = arho(:,:) + tmpPlane(:,:,1)
      au(1,:,:) = au(1,:,:) + tmpPlane(:,:,2)
      au(2,:,:) = au(2,:,:) + tmpPlane(:,:,3)
      au(3,:,:) = au(3,:,:) + tmpPlane(:,:,4)
      ap(:,:)   = ap(:,:)   + tmpPlane(:,:,5)
      aT(:,:)   = aT(:,:)   + tmpPlane(:,:,6)
      amu(:,:)  = amu(:,:)  + tmpPlane(:,:,7)
      aka(:,:)  = aka(:,:)  + tmpPlane(:,:,8)
      aru(1,:,:) = aru(1,:,:) + tmpPlane(:,:,9 ) !ru
      aru(2,:,:) = aru(2,:,:) + tmpPlane(:,:,10 ) !rv
      aru(3,:,:) = aru(3,:,:) + tmpPlane(:,:,11) !rw
      arT(:,:)   = arT(:,:)   + tmpPlane(:,:,12)
      aruu(1,:,:) = aruu(1,:,:) + tmpPlane(:,:,13) !ruu
      aruu(2,:,:) = aruu(2,:,:) + tmpPlane(:,:,14) !ruv
      aruu(3,:,:) = aruu(3,:,:) + tmpPlane(:,:,15) !ruw
      aruu(4,:,:) = aruu(4,:,:) + tmpPlane(:,:,16) !rvv
      aruu(5,:,:) = aruu(5,:,:) + tmpPlane(:,:,17) !rvw
      aruu(6,:,:) = aruu(6,:,:) + tmpPlane(:,:,18) !rww
      aTauij(1,:,:) = aTauij(1,:,:) + tmpPlane(:,:,19) !tau_xx
      aTauij(2,:,:) = aTauij(2,:,:) + tmpPlane(:,:,20) !tau_xy
      aTauij(3,:,:) = aTauij(3,:,:) + tmpPlane(:,:,21) !tau_xz
      aTauij(4,:,:) = aTauij(4,:,:) + tmpPlane(:,:,22) !tau_yy
      aTauij(5,:,:) = aTauij(5,:,:) + tmpPlane(:,:,23) !tau_yz
      aTauij(6,:,:) = aTauij(6,:,:) + tmpPlane(:,:,24) !tau_zz
      aqj(1,:,:) = aqj(1,:,:) + tmpPlane(:,:,25) !q_x
      aqj(2,:,:) = aqj(2,:,:) + tmpPlane(:,:,26) !q_y
      aqj(3,:,:) = aqj(3,:,:) + tmpPlane(:,:,27) !q_z
    enddo
    deallocate(tmpPlane)
  endif
  write(stdout,*)
  write(stdout,*) " averaging done!"
! 1-D and 2-D wall properties
  tauw = aTauij(3,1,:) ! mu*dudy
  qw = aqj(1,1,:) ! absolute value
  rhow = arho(1,:)
  Tw = aT(1,:)
  muw = amu(1,:)*Re
  kaw = aka(1,:)*Re*Pra*Ec
  ! 2-D
  do k=1,xsize(3)  
    do i=1,xsize(1)
        u_tau(i,k) = sqrt(aTauij(3,1,k)/arho(1,k)) ! viscous so only wall units
        u_tau_sl(i,k) = sqrt(aTauij(3,1,k)/arho(i,k)) ! sl hold for semi-local viscous length scale
        Re_tau(i,k) = u_tau(i,k) * arho(1,k)/amu(1,k)  ! to be multiplied by x(2) to obtain y+ 
        Re_tau_sl(i,k) = u_tau_sl(i,k) * arho(i,k)/amu(i,k)  ! to be multiplied by x(2) to obtain y*
    enddo
  enddo
  Re_tauw=Re_tau(1,:)
  ! Momentum thickness
  do k=1,xsize(3)
    call trapzint(xsize(1),x,arho(:,k)*au(3,:,k)*( 1.0_mytype-au(3,:,k) ),theta)
    theta_vector(k)=theta
  enddo
  call assemble_globalz1D(Tw,Tw_global,kstart,kend)
  call assemble_globalz1D(rhow,rhow_global,kstart,kend)
  call assemble_globalz1D(muw,muw_global,kstart,kend)
  call assemble_globalz1D(kaw,kaw_global,kstart,kend)
  call assemble_globalz1D(tauw,tauw_global,kstart,kend)
  call assemble_globalz1D(qw,qw_global,kstart,kend)
  call assemble_globalz1D(Re_tauw,Re_tauw_global,kstart,kend)
  call assemble_globalz1D(theta_vector,theta_global,kstart,kend)
! Perturbation and RMS values
  if (rms_flag==1) then
    write(stdout,*) "o--------------------------------------------------o"
    write(stdout,*) 'RMS'
    write(stdout,*)
    nfiles_rms = (iend_rms-istart_rms)/istep_rms+1
    ! Obtain rms of tauw by reading x-planes
    allocate(tmpPlane (xsize(2),xsize(3),2))
    allocate(tauw_rms(1:partfy%ysz(3)));   tauw_rms = 0.0_mytype;
    allocate(tauw_rms_global (nz_global) )
    ! Check whether the required plane data exists
    write(cha,'(I0.7)') istart_rms
    inquire(file='postproc/planes/xpl.1.mu.'//cha//'.bin', exist=exist1)    ! mu (i=1)
    inquire(file='postproc/planes/xpl.1.strxz.'//cha//'.bin', exist=exist2) ! dudy (i=1)
    if ((exist1) .and. (exist2)) then
      if (index_rms_xpl(1) .eq. 1) then ! planes at wall
        do istep = istart_rms,iend_rms,istep_rms
          write(cha,'(I0.7)') istep
          if (nrank==0) write(stdout,'(A, A10)') 'reading 2D planes for RMS of tauw:    ', cha
            i = index_rms_xpl(1)
            write(cha3,'(I0)') i
            ! Read plane data
            inquire(iolength=lenr) tmpPlane(1,1,1)
            open(10,file='postproc/planes/xpl.'//trim(cha3)//'.mu.'//cha//'.bin',&
                    form='unformatted',status="old",access='direct',recl=xsize(2)*lenr)
            do k=1,xsize(3)
              read(10, rec=k+kstart-1) (tmpPlane(j,k,1), j=1,xsize(2))
            enddo
            close(10)
            inquire(iolength=lenr) tmpPlane(1,1,2)
            open(10,file='postproc/planes/xpl.'//trim(cha3)//'.strxz.'//cha//'.bin',&
                    form='unformatted',status="old",access='direct',recl=xsize(2)*lenr)
            do k=1,xsize(3)
              read(10, rec=k+kstart-1) (tmpPlane(j,k,2), j=1,xsize(2))
            enddo
            close(10)
            do k=1,xsize(3)
              do j=1,xsize(2)
                tauw_local  = tmpPlane(j,k,1)*tmpPlane(j,k,2) ! mu*dudy
                tauw_rms(k) = tauw_rms(k) + (tauw_local - aTauij(3,1,k))*(tauw_local - aTauij(3,1,k))
              enddo
            enddo
        enddo
        tauw_rms = dsqrt(tauw_rms/dble(nfiles_rms*xsize(2)))
      endif
    else
      if (nrank==0) write(*,*) 'RMS of tauw is not obtained'
      tauw_rms = 0.0_mytype
    endif
    deallocate(tmpPlane)
    call assemble_globalz1D(tauw_rms,tauw_rms_global,kstart,kend)
    ! Obtain Quadrant components by reading x-planes
    allocate(tmpPlane (xsize(1),xsize(2),3))
    allocate(aQ4    (5,xsize(1),size(index_rms_zpl))); aQ4     = 0.0_mytype; !1:w'>0 & u'>0 2:w'<0 & u'>0 ...
    allocate(countQ4(4,xsize(1),size(index_rms_zpl))); countQ4 = 0;

    ! Check whether the required plane data exists
    write(cha,'(I0.7)') istart_rms
    write(cha3,'(I0)') index_rms_zpl(1)
    inquire(file='postproc/planes/zpl.'//trim(cha3)//'.u.'//cha//'.bin', exist=exist1)
    inquire(file='postproc/planes/zpl.'//trim(cha3)//'.w.'//cha//'.bin', exist=exist2)
    inquire(file='postproc/planes/zpl.'//trim(cha3)//'.r.'//cha//'.bin', exist=exist3)
    if ((exist1) .and. (exist2) .and. (exist3)) then
      do istep = istart_rms,iend_rms,istep_rms
        write(cha,'(I0.7)') istep
        if (nrank==0) write(stdout,'(A, A10)') 'reading planes for Quadrant analysis: ', cha
        do izpl = 1,size(index_rms_zpl)
          kglobal = index_rms_zpl(izpl)
          if ((kstart .le. kglobal) .and. (kend .ge. kglobal)) then
            k=kglobal-kstart+1
            write(cha3,'(I0)') kglobal
            ! Read plane data
            open(10,file='postproc/planes/zpl.'//trim(cha3)//'.u.'//cha//'.bin',&
                    form='unformatted',status="old",access='stream')
            read(10)tmpPlane(:,:,1)
            close(10)
            open(10,file='postproc/planes/zpl.'//trim(cha3)//'.w.'//cha//'.bin',&
                    form='unformatted',status="old",access='stream')
            read(10)tmpPlane(:,:,2)
            close(10)
            open(10,file='postproc/planes/zpl.'//trim(cha3)//'.r.'//cha//'.bin',&
                    form='unformatted',status="old",access='stream')
            read(10)tmpPlane(:,:,3)
            close(10)
            do i=1,xsize(1)
              do j=1,xsize(2)
                u_flu = tmpPlane(i,j,1)  - aru(1,i,k)/arho(i,k)
                w_flu = tmpPlane(i,j,2)  - aru(3,i,k)/arho(i,k)
                if (w_flu .gt. 0.0_mytype .and. u_flu .gt. 0.0_mytype) then !Q1
                  aQ4(1,i,izpl) = aQ4(1,i,izpl) + tmpPlane(i,j,3)*u_flu*w_flu
                  countQ4(1,i,izpl) = countQ4(1,i,izpl) + 1
                elseif (w_flu .lt. 0.0_mytype .and. u_flu .gt. 0.0_mytype) then !Q2
                  aQ4(2,i,izpl) = aQ4(2,i,izpl) + tmpPlane(i,j,3)*u_flu*w_flu
                  countQ4(2,i,izpl) = countQ4(2,i,izpl) + 1
                elseif (w_flu .lt. 0.0_mytype .and. u_flu .lt. 0.0_mytype) then !Q3
                  aQ4(3,i,izpl) = aQ4(3,i,izpl) + tmpPlane(i,j,3)*u_flu*w_flu
                  countQ4(3,i,izpl) = countQ4(3,i,izpl) + 1
                elseif (w_flu .gt. 0.0_mytype .and. u_flu .lt. 0.0_mytype) then !Q4
                  aQ4(4,i,izpl) = aQ4(4,i,izpl) + tmpPlane(i,j,3)*u_flu*w_flu
                  countQ4(4,i,izpl) = countQ4(4,i,izpl) + 1
                endif
                aQ4(5,i,izpl) = aQ4(5,i,izpl) + tmpPlane(i,j,3)*u_flu*w_flu !Qtotal
              enddo
            enddo
          endif
        enddo
      enddo
      do izpl = 1,size(index_rms_zpl)
        kglobal = index_rms_zpl(izpl)
        if ((kstart .le. kglobal) .and. (kend .ge. kglobal)) then
          do i=1,xsize(1)
            aQ4(1,i,izpl) = aQ4(1,i,izpl)/dble(countQ4(1,i,izpl))
            aQ4(2,i,izpl) = aQ4(2,i,izpl)/dble(countQ4(2,i,izpl))
            aQ4(3,i,izpl) = aQ4(3,i,izpl)/dble(countQ4(3,i,izpl))
            aQ4(4,i,izpl) = aQ4(4,i,izpl)/dble(countQ4(4,i,izpl))
            aQ4(5,i,izpl) = aQ4(5,i,izpl)/dble(xsize(2))/dble(nfiles_rms)
          enddo
        endif
      enddo
    else
      if (nrank==0) write(*,*) 'Quadrant analysis is not conducted'
      aQ4 = 0.0_mytype
    endif
    deallocate(tmpPlane)
  else
      allocate(tauw_rms_global (nz_global) ); tauw_rms_global = 0.0_mytype;
  endif

!===============================================================================================!
!                                         OUTPUTS
!===============================================================================================!
! Store values in simple .txt file
  if (nrank==0) then
    open(18,file = 'postproc/wall_prop.txt')
    write(18,*) 'k - z - Re_delta - Tw - rhow - muw - kaw - Cf - qw - Re_tauw - x^+ - z^+ - Cf_rms'
    k = 1
    write(18,*) k, z_global(k), sqrt((z_global(k)+zStartDNS)*Re), Tw_global(k), rhow_global(k), &
                   muw_global(k), kaw_global(k), &
                   tauw_global(k), qw_global(k), Re_tauw_global(k), Re_tauw_global(k)*x(2), &
                   Re_tauw_global(k)*(z_global(2)-z_global(1)), tauw_rms_global(k)
    do k=2,nz_global
      write(18,*) k, z_global(k), sqrt((z_global(k)+zStartDNS)*Re), Tw_global(k), rhow_global(k), &
                  muw_global(k), kaw_global(k), &
                  tauw_global(k), qw_global(k), Re_tauw_global(k), Re_tauw_global(k)*x(2), &
                  Re_tauw_global(k)*(z_global(k)-z_global(k-1)), tauw_rms_global(k)
    enddo
    close(18)
  endif

! 2D (Wall-normal directions)
  izpl = 1
  kglobal = index_rms_zpl(izpl)
  if ((kstart .le. kglobal) .and. (kend .ge. kglobal)) then
    k=kglobal-kstart+1
    open(18,file = 'postproc/wall_normal_prop.txt')
    write(18,*) 'k index:', kglobal
    write(18,*) 'Re_delta:', sqrt((z_global(kglobal)+zStartDNS)*Re)
    write(18,*) 'i - x - xp - xst - r - u - v - w - p - t - wplus'
    write(18,*) 'ruw/tauw - sqrt(ruu/tauw) - sqrt(rvv/tauw) - sqrt(rww/tauw) - tauxy'
    do i=1,xsize(1)
      write(18,*) i, x(i), Re_tau(i,k)*x(i), Re_tau_sl(i,k)*x(i), arho(i,k), au(1,i,k), &
                  au(2,i,k), au(3,i,k), ap(i,k), aT(i,k), au(3,i,k)/u_tau(1,k), &
                  (aruu(3,i,k) - aru(1,i,k)*aru(3,i,k)/arho(i,k))/aTauij(3,1,k), &
                  dsqrt((aruu(1,i,k) - aru(1,i,k)*aru(1,i,k)/arho(i,k))/aTauij(3,1,k)), &
                  dsqrt((aruu(4,i,k) - aru(2,i,k)*aru(2,i,k)/arho(i,k))/aTauij(3,1,k)), &
                  dsqrt((aruu(6,i,k) - aru(3,i,k)*aru(3,i,k)/arho(i,k))/aTauij(3,1,k)), &
                  aTauij(3,i,k)/aTauij(3,1,k)
    enddo
    close(18)
  endif

! 2D (Mean profiles for recycling/rescaling)
  kglobal = minloc(abs(z_recycle-z_global),1) !Set k-index of your scheduled recycling position (z_recycle)
  if ((kstart .le. kglobal) .and. (kend .ge. kglobal)) then
    k=kglobal-kstart+1
    open(18,file = 'rescale/mean_values.txt')
    write(18,*) 'u,v,w,pre,tem'
    do i=1,xsize(1)
      write(18,*) au(1,i,k),au(2,i,k),au(3,i,k),ap(i,k),aT(i,k)
    enddo
    close(18)
  endif

! Quadrant analysis (wall-normal direction)
  if (rms_flag==1) then
    izpl = 1
    kglobal = index_rms_zpl(izpl)
    if ((kstart .le. kglobal) .and. (kend .ge. kglobal)) then
      k=kglobal-kstart+1
      open(18,file = 'postproc/wall_normal_Q4.txt')
      write(18,*) 'k index:', kglobal
      write(18,*) 'Re_delta:', sqrt((z_global(kglobal)+zStartDNS)*Re)
      write(18,*) 'i - x - xp - xst - ruw/tauw(fav) - Q1 - Q2 - Q3 - Q4 - Qtotal'
      do i=1,xsize(1)
        write(18,*) i, x(i), Re_tau(i,k)*x(i), Re_tau_sl(i,k)*x(i), &
                    (aruu(3,i,k) - aru(1,i,k)*aru(3,i,k)/arho(i,k))/aTauij(3,1,k), &
                    aQ4(1,i,izpl)/aTauij(3,1,k), aQ4(2,i,izpl)/aTauij(3,1,k), &
                    aQ4(3,i,izpl)/aTauij(3,1,k), aQ4(4,i,izpl)/aTauij(3,1,k), &
                    aQ4(5,i,izpl)/aTauij(3,1,k)
      enddo
    endif
  endif
  if (nrank==0) then
    write(stdout,*)
    write(stdout,*) "o--------------------------------------------------o"
    write(stdout,*)
    write(stdout,*) 'Output properties done!'
  endif

! Spanwise and time FFT
  if (fft_flag==1)  then 
    if (nrank == 0) then
      write(stdout,*) "o--------------------------------------------------o"
      write(stdout,*) 'FFT'
      write(stdout,*)
      nfiles_fft = (iend_fft-istart_fft)/istep_fft+1
      write(stdout,'(A, I10)') 'number of files for FFT:              ', nfiles_fft-1
    endif 
    if (avg_flag == 'restart') then
      do ii_index=1,count
        do k=1,xsize(3)
          do i=1,xsize(1)
            !arho_fluc(ii_index,i,1:xsize(2),k)=arho_time(ii_index,i,1:xsize(2),k)-arho(i,k)
            !au_fluc(ii_index,i,1:xsize(2),k)=au_time(ii_index,i,1:xsize(2),k)-au(1,i,k)
            !av_fluc(ii,i,1:xsize(2),k)=av_time(ii,i,1:xsize(2),k)-au(2,i,k)
            aw_fluc_FFT(ii_index,i,1:xsize(2),k)=aw_time(ii_index,i,1:xsize(2),k)-w_bl(i,1:xsize(2),k)
            !aT_fluc(ii,i,1:xsize(2),k)=aT_time(ii,i,1:xsize(2),k)-aT(i,k)
            aP_fluc_FFT(ii_index,i,1:xsize(2),k)=aP_time(ii_index,i,1:xsize(2),k)-pre_bl(i,1:xsize(2),k)
            arw_fluc_FFT(ii_index,i,1:xsize(2),k)=arw_time(ii_index,i,1:xsize(2),k) - &
                                                  (rho_bl(i,1,k) * w_bl(i,1,k))
          enddo
        enddo
        write(cha,'(I0.7)') istart_pp+istep_pp*(ii_index - 1)
        if (nrank==0) write(*,*) cha
        ! call the FFT routine
        !call spectray(2,spec_rho,jmax,part1,partfy,arho_fluc(ii_index,1:xsize(1),1:xsize(2),1:xsize(3)))
        call spectray(2,spec_w,jmax,part1,partfy,aw_fluc_FFT(ii_index,1:xsize(1),1:xsize(2),1:xsize(3))) 
        !call spectray(2,spec_u,jmax,part1,partfy,au_fluc(ii_index,1:xsize(1),1:xsize(2),1:xsize(3)))
        call spectray(2,spec_rw,jmax,part1,partfy,arw_fluc_FFT(ii_index,1:xsize(1),1:xsize(2),1:xsize(3)))
        call spectray(2,spec_p,jmax,part1,partfy,aP_fluc_FFT(ii_index,1:xsize(1),1:xsize(2),1:xsize(3))) 

    !===============================================================================================!
    !                                         OUTPUT
    !===============================================================================================!
        ! Write out 2D-planes of span and time FFT quantities (written in /postproc/results/fft)
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
    else if (avg_flag == 'planes') then
      allocate(tmpPlane (xsize(2),xsize(3),3))
      allocate(arho_fluc(xsize(1),xsize(2),xsize(3)));  arho_fluc = 0.0_mytype;
      allocate(aw_fluc  (xsize(1),xsize(2),xsize(3)));  aw_fluc = 0.0_mytype;
      allocate(au_fluc  (xsize(1),xsize(2),xsize(3)));  au_fluc = 0.0_mytype;
      if (nrank==0) write(stdout,'(A, I10)') 'fft_xplane:                           ', index_fft_xpl
      write(cha3,'(I0)') index_fft_xpl
      count=0
      do istep = istart_fft, (iend_fft-istep_fft), istep_fft
        count=count+1
        write(cha,'(I0.7)') istep
        if (nrank==0) write(stdout,'(A, A10)') 'reading 2D planes for FFT:            ', cha
        ! Read plane data
        inquire(iolength=lenr) tmpPlane(1,1,1)
        open(10,file='postproc/planes/xpl.'//trim(cha3)//'.r.'//cha//'.bin',&
                form='unformatted',status="old",access='direct',recl=xsize(2)*lenr)
        do k=1,xsize(3)
          read(10, rec=k+kstart-1) (tmpPlane(j,k,1), j=1,xsize(2))
        enddo
        close(10)
        inquire(iolength=lenr) tmpPlane(1,1,2)
        open(10,file='postproc/planes/xpl.'//trim(cha3)//'.w.'//cha//'.bin',&
                form='unformatted',status="old",access='direct',recl=xsize(2)*lenr)
        do k=1,xsize(3)
          read(10, rec=k+kstart-1) (tmpPlane(j,k,2), j=1,xsize(2))
        enddo
        close(10)
        inquire(iolength=lenr) tmpPlane(1,1,3)
        open(10,file='postproc/planes/xpl.'//trim(cha3)//'.u.'//cha//'.bin',&
                form='unformatted',status="old",access='direct',recl=xsize(2)*lenr)
        do k=1,xsize(3)
          read(10, rec=k+kstart-1) (tmpPlane(j,k,3), j=1,xsize(2))
        enddo
        close(10)
        ! Calculate fluctuations
        do k=1,xsize(3)
          do i=1,xsize(1)
            ! At index_fft_xpl in wall-normal direction 
            arho_fluc(i,1:xsize(2),k) =tmpPlane(1:xsize(2),k,1)-rho_bl(index_fft_xpl,1:xsize(2),k)!arho(index_fft_xpl,k)
            aw_fluc(i,1:xsize(2),k)   =tmpPlane(1:xsize(2),k,2)-w_bl(index_fft_xpl,1:xsize(2),k)
            au_fluc(i,1:xsize(2),k)   =tmpPlane(1:xsize(2),k,3)-u_bl(index_fft_xpl,1:xsize(2),k)
          enddo
        enddo
        ! FFT in spanwise direction
        call spectray(2,spec_rho,jmax,part1,partfy,arho_fluc(1:xsize(1),1:xsize(2),1:xsize(3)))
        call spectray(2,spec_w,jmax,part1,partfy,aw_fluc(1:xsize(1),1:xsize(2),1:xsize(3))) 
        call spectray(2,spec_u,jmax,part1,partfy,au_fluc(1:xsize(1),1:xsize(2),1:xsize(3)))
        do j=1,size_spany
          write(cha2,'(I0)') index_fft_span(j)     
          tmp_complex = spec_rho(1:xsize(1),1:xsize(2),1:xsize(3))
          call decomp_2d_write_plane(1,tmp_complex,2,index_fft_span(j)+1, &
                                     '.','postproc/results/fft/fft_r_span'//cha2//'_'//cha//'.bin', 'dummy')
          tmp_complex = spec_w(1:xsize(1),1:xsize(2),1:xsize(3))
          call decomp_2d_write_plane(1,tmp_complex,2,index_fft_span(j)+1, &
                                     '.','postproc/results/fft/fft_w_span'//cha2//'_'//cha//'.bin', 'dummy')
          tmp_complex = spec_u(1:xsize(1),1:xsize(2),1:xsize(3))
          call decomp_2d_write_plane(1,tmp_complex,2,index_fft_span(j)+1, &
                                     '.','postproc/results/fft/fft_u_span'//cha2//'_'//cha//'.bin', 'dummy')
        enddo
      enddo
      deallocate(tmpPlane)
    end if
    if (nrank==0) write(stdout,*)
    if (nrank==0) write(stdout,*) 'FFT done!'
  endif
  if (nrank == 0) write(stdout,*) "o--------------------------------------------------o"
  if (nrank == 0) write(stdout,*) 'POST done!'
! print total time
  if (nrank == 0) print '("Total time = ",f10.3," minutes.")', (MPI_WTIME() - wt_start)/60.0

!===============================================================================================!
!                                       DEALLOCATION
!===============================================================================================!
    deallocate(rho)
    deallocate(u)
    deallocate(v)
    deallocate(w)
    deallocate(ien)
    deallocate(pre)
    deallocate(tem)
    deallocate(mu)
    deallocate(ka)
  if (avg_flag == 'restart') then
    deallocate(Cp)
    deallocate(qvort)
    deallocate(dilla2)
    deallocate(sxx)
    deallocate(sxy)
    deallocate(sxz)
    deallocate(syy)
    deallocate(syz)
    deallocate(szz)
  endif
  if (perBC(3) .eqv. .false.) then
    deallocate(rho_bl)
    deallocate(u_bl)
    deallocate(v_bl)
    deallocate(w_bl)
    deallocate(ien_bl)
    deallocate(pre_bl)
    deallocate(tem_bl)
    deallocate(mu_bl)
    deallocate(ka_bl)
  endif
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call mpi_finalize(ierr)
end program

! assemble global variables in streamwise direction
  subroutine assemble_globalz1D(local1D,global1D,kstart,kend)
    use decomp_2d
    use mpi 
    implicit none
    integer, intent(IN) :: kstart,kend
    real(mytype), dimension(:), intent(IN) :: local1D 
    real(mytype), dimension(nz_global), intent(OUT) :: global1D 
    real(mytype), allocatable, dimension(:) :: rbuf_1D 
    integer, dimension(3) :: sbuf1D, rbuf1D
    integer :: ierror,k,m,k1,k2, count
    integer, dimension(MPI_STATUS_SIZE) :: status
    if (nrank==0) then
       ! master writes its own data to a global array
          k1 = kstart
          k2 = kend
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
          sbuf1D(1) = kstart
          sbuf1D(2) = kend
          sbuf1D(3) = kend
          count = kend-kstart+1
       ! send partition information
       CALL MPI_SEND(sbuf1D,3,MPI_INTEGER,0,nrank,MPI_COMM_WORLD,ierror)
       ! send data array
       CALL MPI_SEND(local1D,count,real_type,0, &
            nrank+nproc,MPI_COMM_WORLD,ierror) 
    end if
    return
  end subroutine assemble_globalz1D
