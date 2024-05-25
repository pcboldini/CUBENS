
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
  integer :: kstart,kend,klocal,kglobal,iname,lenr
  integer :: istep,nfiles_rms,ixpl,izpl,nfiles_fft,count,size_spany
  TYPE (DECOMP_INFO) :: part1,partfx,partfz,partfy,partftime !!!
  TYPE (DECOMP_INFO) :: partinterp
  character*7 :: cha
  character*1 :: cha2
  character(len=1024) :: cha3
  integer, dimension(3) :: fft_start, fft_end, fft_size
  character(len=30) :: date
  real(8) :: wt_start
  logical :: exist1,exist2,exist3


  ! SCRINS Version number
   real(mytype), parameter                    :: version = 1.1                      ! SCRINS version

  ! underlying FFT library only needs to be initialised once
  !logical, save :: fft_initialised = .false.

  real(mytype) :: factAvg

  !read & write
  real(mytype), allocatable, dimension(:,:,:)    :: tmpPlane           !(i, j, k) or (i, k, number of planes)
  complex(mytype), allocatable, dimension(:,:,:) :: tmp_complex        !(i, j, k)

  !span_ave (stats)
  character(5) :: stat_name(26)
  real(mytype), allocatable, dimension(:,:)   ::  arho, ap, aT, amu, arT
  real(mytype), allocatable, dimension(:,:,:) ::  au, aru, aruu, aTauij, aqj

  real(mytype), allocatable, dimension(:) :: Tw, rhow, tauw, qw, Re_tauw
  real(mytype), allocatable, dimension(:) :: Tw_global, rhow_global, tauw_global, qw_global, Re_tauw_global

  real(mytype), allocatable, dimension(:,:) :: u_tau, u_tau_sl
  real(mytype), allocatable, dimension(:,:) :: Re_tau, Re_tau_sl

  !RMS (tauw)
  real(mytype), allocatable, dimension(:)   ::  tauw_rms, tauw_rms_global
  real(mytype) :: tauw_local

  !RMS (quadrant)
  real(mytype), allocatable, dimension(:,:,:) :: aQ4
  integer, allocatable, dimension(:,:,:) :: countQ4
  real(mytype) :: u_flu, w_flu

  !fft
  real(mytype), allocatable, dimension(:,:,:) :: arho_fluc
  real(mytype), allocatable, dimension(:,:,:) :: aw_fluc
  real(mytype), allocatable, dimension(:,:,:) :: au_fluc
  complex(mytype), allocatable, dimension(:,:,:) :: spec_rho
  complex(mytype), allocatable, dimension(:,:,:) :: spec_w
  complex(mytype), allocatable, dimension(:,:,:) :: spec_u

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: kk = 6
  real(mytype) :: cff
  real(mytype), dimension(kk) :: in1, out1
  real(mytype), dimension(kk) :: test_in, test_out, test_out1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

interface
   subroutine assemble_globalz1D(local1D,global1D,kstart,kend)
     use decomp_2d
     integer, intent(IN) :: kstart,kend
     real(mytype), dimension(:), intent(IN) :: local1D !!!!
     real(mytype), dimension(nz_global), intent(OUT) :: global1D !!!!
   end subroutine assemble_globalz1D
end interface


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
  select case (INIT)
    case("initBoundaryLayer")
      call read_initBL_params()
  end select

  nfiles_rms = (iend_rms-istart_rms)/istep_rms+1
  nfiles_fft = (iend_fft-istart_fft)/istep_fft+1

! Init MPI and decomp_2d
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

 call decomp_info_init(nfiles_fft/2+1,imax,kmax,partftime)
 call decomp_info_init(imax,jmax,kmax/2+1,partfz)
 call decomp_info_init(imax,jmax/2+1,kmax,partfy)
 call decomp_info_init(imax/2+1,jmax,kmax,partfx)

if (nrank==0) write(stdout,*) 'xpencil (partftime%xsz)', partfy%xsz
if (nrank==0) write(stdout,*) 'ypencil (partftime%ysz)', partfy%ysz
if (nrank==0) write(stdout,*) 'zpencil (partftime%zsz)', partfy%zsz

if (nrank==0) write(stdout,*) 'stats_step', stats_step
if (nrank==0) write(stdout,*) 'stats_time_rate', stats_time_rate

if (rms_flag==1) then
  if (nrank==0) write(stdout,*) 'istart_rms', istart_rms
  if (nrank==0) write(stdout,*) 'iend_rms', iend_rms
  if (nrank==0) write(stdout,*) 'index_rms_xpl', index_rms_xpl
  if (nrank==0) write(stdout,*) 'index_rms_zpl', index_rms_zpl
endif

if (fft_flag==1) then 
  if (nrank==0) write(stdout,*) 'index_fft_span', index_fft_span
  if (nrank==0) write(stdout,*) 'fft_xplane', index_fft_xpl
  size_spany=size(index_fft_span)
else
  index_fft_span=1
  size_spany=1
endif

! allocate main variables
  allocate(arho(  xsize(1),xsize(3)));     arho = 0.0_mytype;
  allocate(  ap(  xsize(1),xsize(3)));     ap   = 0.0_mytype;
  allocate(  aT(  xsize(1),xsize(3)));     aT   = 0.0_mytype;
  allocate( amu(  xsize(1),xsize(3)));     amu  = 0.0_mytype;
  allocate( arT(  xsize(1),xsize(3)));     arT  = 0.0_mytype;

  allocate(    au(3,xsize(1),xsize(3)));   au     = 0.0_mytype;  !3 components
  allocate(   aru(3,xsize(1),xsize(3)));   aru    = 0.0_mytype;  !3 components
  allocate(  aruu(6,xsize(1),xsize(3)));   aruu   = 0.0_mytype;  !6 components
  allocate(aTauij(6,xsize(1),xsize(3)));   aTauij = 0.0_mytype;  !6 components
  allocate(   aqj(3,xsize(1),xsize(3)));   aqj    = 0.0_mytype;  !3 components
   
  allocate(Tw     (1:partfy%ysz(3))    )
  allocate(rhow   (1:partfy%ysz(3))    )
  allocate(tauw   (1:partfy%ysz(3))    )
  allocate(qw     (1:partfy%ysz(3))    )
  allocate(Re_tauw(1:partfy%ysz(3))    )

  allocate(Tw_global     (nz_global)   )
  allocate(rhow_global   (nz_global)   )
  allocate(tauw_global   (nz_global)   )
  allocate(qw_global     (nz_global)   )
  allocate(Re_tauw_global(nz_global)   )

  allocate(u_tau    (xsize(1),xsize(3)));   u_tau     = 0.0_mytype;
  allocate(u_tau_sl (xsize(1),xsize(3)));   u_tau_sl  = 0.0_mytype;
  allocate(Re_tau   (xsize(1),xsize(3)));   Re_tau    = 0.0_mytype;
  allocate(Re_tau_sl(xsize(1),xsize(3)));   Re_tau_sl = 0.0_mytype;

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
  call init_PARAM()

  if (nrank==0) write(stdout,*) 'initializing grid'
  call init_grid()

   if (nrank==0) write(stdout,*) 'initializing finite difference coeffs'
  call init_derivCoeffs()

  if (nrank==0) write(stdout,*) 'initializing rhs'
  call init_rhs()

  if (nrank==0) write(stdout,*) 'p_row_pp', p_row_pp
  if (nrank==0) write(stdout,*) 'p_col_pp', p_col_pp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------
! SPANWISE AVERAGING 
!----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  count=0
  wt_start = MPI_WTIME()

  !Check number of stats
  if (nrank==0) write(*,*) 'number of stats files: ', size(stats_step)
  if (size(stats_step) .ne. size(stats_time_rate)) then
     if(nrank==0) write(stdout,*) 'Check number of stats in config.h!'
     call MPI_FILE_CLOSE(fh,ierr)
     call decomp_2d_finalize
     call MPI_FINALIZE(ierr)
     stop
  end if
  factAvg = 1.0_mytype / sum(stats_time_rate)

  !Get partition information
  do k=1,nz_global
     if (z(1) .eq. z_global(k)) kstart = k
     if (z(xsize(3)) .eq. z_global(k)) kend = k
  enddo

  !Set file name (number of elements are decided when the variable is defined)
  stat_name = ['r    ','u    ','v    ','w    ','p    ','t    ','mu   ', &
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
     if (nrank==0) write(*,'(/,A,A)') 'reading stats: ', cha

     do iname = 1,size(stat_name)
        inquire(iolength=lenr) tmpPlane(1,1,iname)
        open(10,file='postproc/planes/stats.'//trim(stat_name(iname))//'.'//cha//'.bin',&
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

     aru(1,:,:) = aru(1,:,:) + tmpPlane(:,:,8 ) !ru
     aru(2,:,:) = aru(2,:,:) + tmpPlane(:,:,9 ) !rv
     aru(3,:,:) = aru(3,:,:) + tmpPlane(:,:,10) !rw
     arT(:,:)   = arT(:,:)   + tmpPlane(:,:,11)

     aruu(1,:,:) = aruu(1,:,:) + tmpPlane(:,:,12) !ruu
     aruu(2,:,:) = aruu(2,:,:) + tmpPlane(:,:,13) !ruv
     aruu(3,:,:) = aruu(3,:,:) + tmpPlane(:,:,14) !ruw
     aruu(4,:,:) = aruu(4,:,:) + tmpPlane(:,:,15) !rvv
     aruu(5,:,:) = aruu(5,:,:) + tmpPlane(:,:,16) !rvw
     aruu(6,:,:) = aruu(6,:,:) + tmpPlane(:,:,17) !rww

     aTauij(1,:,:) = aTauij(1,:,:) + tmpPlane(:,:,18) !tau_xx
     aTauij(2,:,:) = aTauij(2,:,:) + tmpPlane(:,:,19) !tau_xy
     aTauij(3,:,:) = aTauij(3,:,:) + tmpPlane(:,:,20) !tau_xz
     aTauij(4,:,:) = aTauij(4,:,:) + tmpPlane(:,:,21) !tau_yy
     aTauij(5,:,:) = aTauij(5,:,:) + tmpPlane(:,:,22) !tau_yz
     aTauij(6,:,:) = aTauij(6,:,:) + tmpPlane(:,:,23) !tau_zz

     aqj(1,:,:) = aqj(1,:,:) + tmpPlane(:,:,24) !q_x
     aqj(2,:,:) = aqj(2,:,:) + tmpPlane(:,:,25) !q_y
     aqj(3,:,:) = aqj(3,:,:) + tmpPlane(:,:,26) !q_z

  enddo
  deallocate(tmpPlane)

  !Obtain wall quantities
  Tw   = aT(1,:)
  rhow = arho(1,:)
  tauw = aTauij(3,1,:) ! mu*dudy
  qw   = aqj(1,1,:) ! absolute value

  do k=1,xsize(3)
   do i=1,xsize(1)
       u_tau(i,k) = sqrt(aTauij(3,1,k)/arho(1,k))    ! viscous so only wall units
       u_tau_sl(i,k) = sqrt(aTauij(3,1,k)/arho(i,k)) ! sl hold for semi-local viscous length scale

       Re_tau(i,k) = u_tau(i,k) * arho(1,k)/amu(1,k)        !! to be multiplied by x(2) to obtain y+ 
       Re_tau_sl(i,k) = u_tau_sl(i,k) * arho(i,k)/amu(i,k)  !! to be multiplied by x(2) to obtain y*
   enddo
  enddo

  Re_tauw=Re_tau(1,:)

  call assemble_globalz1D(Tw,Tw_global,kstart,kend)
  call assemble_globalz1D(rhow,rhow_global,kstart,kend)
  call assemble_globalz1D(tauw,tauw_global,kstart,kend)
  call assemble_globalz1D(qw,qw_global,kstart,kend)
  call assemble_globalz1D(Re_tauw,Re_tauw_global,kstart,kend)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------
! RMS & Quadrant analysis
!----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rms_flag==1) then
     !---------------------------------------
     !Obtain rms of tauw by reading x-planes
     !---------------------------------------
     allocate(tmpPlane (xsize(2),xsize(3),2))
     allocate(tauw_rms(1:partfy%ysz(3)));   tauw_rms = 0.0_mytype;
     allocate(tauw_rms_global (nz_global) )

     !Check whether the required plane data exists
     write(cha,'(I0.7)') istart_rms
     inquire(file='postproc/planes/xpl.1.mu.'//cha//'.bin', exist=exist1)    !mu (i=1)
     inquire(file='postproc/planes/xpl.1.strxz.'//cha//'.bin', exist=exist2) !dudy (i=1)
     if ((exist1) .and. (exist2)) then
     
        if (index_rms_xpl(1) .eq. 1)then !planes at wall
           do istep = istart_rms,iend_rms,istep_rms
              write(cha,'(I0.7)') istep
              if (nrank==0) write(*,'(/,A,A)') 'reading 2D planes for RMS of tauw: ', cha

              i = index_rms_xpl(1)
              write(cha3,'(I0)') i

              !Read plane data------
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
              !---------------------

              do k=1,xsize(3)
                 do j=1,xsize(2)
                    tauw_local  = tmpPlane(j,k,1)*tmpPlane(j,k,2) !mu*dudy
                    tauw_rms(k) = tauw_rms(k) + (tauw_local - aTauij(3,1,k))*(tauw_local - aTauij(3,1,k))
                 enddo
              enddo
           enddo
           tauw_rms = dsqrt(tauw_rms/dble(nfiles_rms*xsize(2)))
        endif
     else
        if (nrank==0) write(*,*) 'RMS of tauw did not be obtained'
        tauw_rms = 0.0_mytype
     endif
     deallocate(tmpPlane)
     call assemble_globalz1D(tauw_rms,tauw_rms_global,kstart,kend)

     !-----------------------------------------------
     !Obtain Quadrant components by reading x-planes
     !-----------------------------------------------
     allocate(tmpPlane (xsize(1),xsize(2),3))
     allocate(aQ4    (5,xsize(1),size(index_rms_zpl))); aQ4     = 0.0_mytype; !1:w'>0 & u'>0 2:w'<0 & u'>0 ...
     allocate(countQ4(4,xsize(1),size(index_rms_zpl))); countQ4 = 0;

     !Check whether the required plane data exists
     write(cha,'(I0.7)') istart_rms
     write(cha3,'(I0)') index_rms_zpl(1)
     inquire(file='postproc/planes/zpl.'//trim(cha3)//'.u.'//cha//'.bin', exist=exist1)
     inquire(file='postproc/planes/zpl.'//trim(cha3)//'.w.'//cha//'.bin', exist=exist2)
     inquire(file='postproc/planes/zpl.'//trim(cha3)//'.r.'//cha//'.bin', exist=exist3)
     if ((exist1) .and. (exist2) .and. (exist3)) then

        do istep = istart_rms,iend_rms,istep_rms
           write(cha,'(I0.7)') istep
           if (nrank==0) write(*,'(/,A,A)') 'reading 2D planes for Quadrant analysis: ', cha

           do izpl = 1,size(index_rms_zpl)
              kglobal = index_rms_zpl(izpl)
              if ((kstart .le. kglobal) .and. (kend .ge. kglobal)) then
                 k=kglobal-kstart+1
                 write(cha3,'(I0)') kglobal

                 !Read plane data------
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
                 !---------------------

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
        if (nrank==0) write(*,*) 'Quadrant analysis did not be conducted'
        aQ4 = 0.0_mytype
     endif
     deallocate(tmpPlane)
  else
     allocate(tauw_rms_global (nz_global) ); tauw_rms_global = 0.0_mytype;
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------
! Output data
!----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!! 2D (Wall planes)
  if (nrank==0) then
     open(18,file = 'postproc/wall_prop.txt')
     write(18,*) 'k - z - Re_delta - Tw - rhow - Cf - qw - Re_tauw - x^+ - z^+ - Cf_rms'
     k = 1
     write(18,*) k, z_global(k), sqrt((z_global(k)+zStartDNS)*Re), Tw_global(k), rhow_global(k), &
                 tauw_global(k), qw_global(k), Re_tauw_global(k), Re_tauw_global(k)*x(2), &
                 Re_tauw_global(k)*(z_global(2)-z_global(1)), tauw_rms_global(k)
     do k=2,nz_global
        write(18,*) k, z_global(k), sqrt((z_global(k)+zStartDNS)*Re), Tw_global(k), rhow_global(k), &
                    tauw_global(k), qw_global(k), Re_tauw_global(k), Re_tauw_global(k)*x(2), &
                    Re_tauw_global(k)*(z_global(k)-z_global(k-1)), tauw_rms_global(k)
     enddo
     close(18)
  endif

  !!!! 2D (Wall-normal directions)
  kglobal = 1668 !index_rms_zpl(1)
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

  !!!! 2D (Mean profiles for RR)
  kglobal = 1668 !Set k-index of your scheduled recycling position
  if ((kstart .le. kglobal) .and. (kend .ge. kglobal)) then
     k=kglobal-kstart+1
     open(18,file = 'postproc/Mean_RR.txt')
     write(18,*) 'u,v,w,pre,tem'
     do i=1,xsize(1)
        write(18,*) au(1,i,k),au(2,i,k),au(3,i,k),ap(i,k),aT(i,k)
     enddo
     close(18)
  endif

  !!!! Quadrant analysis (Wall-normal directions)
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

     if (nrank==0) write(*,*) 'number of files for FFT: ', nfiles_fft-1

     !Allocation
     allocate(tmpPlane (xsize(2),xsize(3),3))
     allocate( arho_fluc(xsize(1),xsize(2),xsize(3)));  arho_fluc = 0.0_mytype;
     allocate( aw_fluc  (xsize(1),xsize(2),xsize(3)));  aw_fluc = 0.0_mytype;
     allocate( au_fluc  (xsize(1),xsize(2),xsize(3)));  au_fluc = 0.0_mytype;
     allocate( spec_rho (1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3)))
     allocate( spec_w   (1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3)))
     allocate( spec_u   (1:partfy%ysz(1), 1:partfy%ysz(2), 1:partfy%ysz(3)))
     allocate( tmp_complex (xsize(1),xsize(2),xsize(3)))

     write(cha3,'(I0)') index_fft_xpl

     count=0
     do istep = istart_fft, (iend_fft-istep_fft), istep_fft
        count=count+1
        write(cha,'(I0.7)') istep
        if (nrank==0) write(*,'(/,A,A)') 'reading 2D planes for FFT: ', cha

        !Read plane data------
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
        !---------------------

        !Calc fluctuation
        do k=1,xsize(3)
           do i=1,xsize(1)
              !A same value in wall-normal direction
              arho_fluc(i,1:xsize(2),k)=tmpPlane(1:xsize(2),k,1)-arho(index_fft_xpl,k)
              aw_fluc(i,1:xsize(2),k)  =tmpPlane(1:xsize(2),k,2)-au(3,index_fft_xpl,k)
              au_fluc(i,1:xsize(2),k)  =tmpPlane(1:xsize(2),k,3)-au(1,index_fft_xpl,k)
           enddo
        enddo

        !FFT in spanwise direction
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

     if (nrank==0) write(*,*) 'FFT done!'
     deallocate(tmpPlane)

  end if

  if (nrank == 0) write(*,*)
  if (nrank == 0) write (stdout, *) "o----------------------------------------------------------------------------------o"
  if (nrank == 0) write(*,*) 'POST done!'
  if (nrank == 0) print '("Total time = ",f10.3," minutes.")', (MPI_WTIME() - wt_start)/60.0

  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call mpi_finalize(ierr)

end program


subroutine assemble_globalz1D(local1D,global1D,kstart,kend)
  
  use decomp_2d
  use mpi
  
  implicit none
  
  integer, intent(IN) :: kstart,kend
  real(mytype), dimension(:), intent(IN) :: local1D !!!!
  real(mytype), dimension(nz_global), intent(OUT) :: global1D !!!!
  
  real(mytype), allocatable, dimension(:) :: rbuf_1D !!!!!
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
          nrank+nproc,MPI_COMM_WORLD,ierror) !!!!!
  end if
  
  return
end subroutine assemble_globalz1D

