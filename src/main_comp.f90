! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! main module

program DNS
  use decomp_2d
  use decomp_2d_io
  use io_std_units
  use mpi
  use mod_param
  use mod_halo
  use mod_grid
  use mod_init
  use mod_eos
  use mod_eos_var
  use mod_eos_visc
  use mod_solve
  use mod_finitediff
  use mod_rhs
  use mod_auxl
  use mod_perturbation
  use mod_boundary
  use iso_fortran_env
  use mod_math
  !@acc use openacc 
  implicit none


!===============================================================================================!
!                                       CUBENS main data
!===============================================================================================!
! Non-conservative, time variables                          
  integer :: ierr,k,j,i
  integer :: istep,istart,countAvg
  real(8) :: wt_start, wt_tmp 
  real(mytype) :: dt, time, velwb, CFL_new, factAvg
  real(mytype), allocatable, dimension(:,:,:) :: rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz,strxz
  real(mytype), allocatable, dimension(:,:,:) :: rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl
#if defined(BL) || defined(TGV) || defined(DHC)
  ! for boundary layer, TGV, and DHC
  real(mytype), allocatable, dimension(:,:,:) :: qave,qtime
#elif defined(CHA)
  real(mytype), allocatable, dimension(:,:) :: qave,qtime
#endif
  TYPE (DECOMP_INFO) :: part
  CHARACTER(100) :: commLinePar1
  CHARACTER(100) :: commLinePar2
  character(len=30) :: date
  character*4 :: cha2
! CUBENS Version number
  real(mytype), parameter                    :: version = 1.0     


!===============================================================================================!
!                                       INITIALIZATION
!===============================================================================================!
! generate today's date
  call fdate(date) 

!---------------------------------------!
!  Read the config.h and init_BL files  !
!---------------------------------------!
! Reading parameters
  call read_config()
#if defined(BL) || defined(CHA)
  call read_init_params()
  ! Writing parameters for user variation 
  call io_writeParams_variation()
#endif

!-----------------------------!
!         MPI library         !
!-----------------------------!
! users can specify the number of cores (p_row and p_col) also from the command line
  if (COMMAND_ARGUMENT_COUNT().eq.2) then
    call GET_COMMAND_ARGUMENT(1, commLinePar1)
    call GET_COMMAND_ARGUMENT(2, commLinePar2)
    read(commLinePar1,*) p_row
    read(commLinePar2,*) p_col
  endif
! init MPI and decomp_2d
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,nrank,ierr)
  call decomp_2d_init(imax,jmax,kmax,p_row,p_col)
  call get_decomp_info(part)
  call comm_init(nrank,p_row,p_col,DECOMP_2D_COMM_CART_X,xsize,imax,jmax,kmax)
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
    write (stdout, '(A,I0,A,I0,A,I0,A)') " Using ", nproc, " MPI processes, ", p_col, " number of procs in z and ", p_row, " in y " 
    write (stdout, *)
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *) "|                     AM MAIN                      |"
    write (stdout, *) "o--------------------------------------------------o"
    write (stdout, *)
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
! initialize Runge-Kutta time integration
  if (nrank==0) then
    write(stdout,*) 'Initializing RK3'
    write(stdout,*)
  endif
  call init_rk3()
! initialize perturbation for transitional boundary layer
  if (nrank==0 .and. pert_calc==1) then
    write(stdout,*) 'Initializing perturbation'
    write(stdout,*)
  endif
  call init_pert()
! initialize boundary conditions
  if (nrank==0) then
    write(stdout,*) 'Initializing boundary conditions'
    write(stdout,*)
  endif
  call init_BC() 
  ! Set communicator for rescale-reintroduction
  if (BC_inl_rescale .eqv. .true.) then
    call comm_init_rescale(z,xsize)
    if (nrank==0) then
      write(stdout,*) 'Initialize communication for recycling and rescaling'
      write(stdout,*)
    endif
  endif


!===============================================================================================!
!                                       ALLOCATION
!===============================================================================================!
! allocate main variables 
  allocate(rho  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(u    (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(v    (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(w    (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(ien  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(pre  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(tem  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(mu   (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(ka   (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(vortx(1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(vorty(1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(vortz(1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  allocate(strxz(1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
  ! save 51 statistics
#if defined(BL) || defined(TGV) || defined(DHC)
  ! for boundary layer, TGV, and DHC
  allocate(qave(xsize(1), xsize(3), 51))
  allocate(qtime(xsize(2), xsize(3), 1))
#elif defined(CHA)
  allocate(qave(xsize(1), 51))
  allocate(qtime(1, 1))
#endif
  rho      = 1.0e30_mytype
  u        = 1.0e30_mytype
  v        = 1.0e30_mytype
  w        = 1.0e30_mytype
  ien      = 1.0e30_mytype
  pre      = 1.0e30_mytype
  tem      = 1.0e30_mytype
  mu       = 1.0e30_mytype
  ka       = 1.0e30_mytype
  vortx    = 1.0e30_mytype
  vorty    = 1.0e30_mytype
  vortz    = 1.0e30_mytype
  strxz    = 1.0e30_mytype 
  qave     = 0.0_mytype
  qtime    = 0.0_mytype
  factAvg  = 0.0_mytype
  countAvg = 0
! copy common arrays to device
!$acc enter data copyin(rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz,strxz,qave,qtime)

! calculate sponge
#if defined(BL)
  if (perBC(3) .eqv. .false.) then
    allocate(rho_bl  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(u_bl    (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(v_bl    (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(w_bl    (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(ien_bl  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(pre_bl  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(tem_bl  (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(mu_bl   (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(ka_bl   (1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
      call initField_BL(part,x,z,rho,u,v,w,ien,pre,tem,mu,ka)
  !$acc update host(rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz,strxz)  
      r_ref(1:xsize(1),1:xsize(3)) = rho(1:xsize(1),1,1:xsize(3))
     ru_ref(1:xsize(1),1:xsize(3)) = rho(1:xsize(1),1,1:xsize(3))*u(1:xsize(1),1,1:xsize(3))
     rv_ref(1:xsize(1),1:xsize(3)) = rho(1:xsize(1),1,1:xsize(3))*v(1:xsize(1),1,1:xsize(3))
     rw_ref(1:xsize(1),1:xsize(3)) = rho(1:xsize(1),1,1:xsize(3))*w(1:xsize(1),1,1:xsize(3))
    ret_ref(1:xsize(1),1:xsize(3)) = rho(1:xsize(1),1,1:xsize(3))*( ien(1:xsize(1),1,1:xsize(3)) &
           + 0.5*(u(1:xsize(1),1,1:xsize(3))**2 + v(1:xsize(1),1,1:xsize(3))**2 + w(1:xsize(1),1,1:xsize(3))**2) )
    p_ref(1:xsize(1),1:xsize(3)) = pre(1:xsize(1),1,1:xsize(3))
#if defined(BL)
    !$acc enter data copyin(rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl)
    !$acc enter data copyin(r_ref,ru_ref,rv_ref,rw_ref,ret_ref,p_ref)
#endif
  endif
#endif 


!===============================================================================================!
!                                       INITIAL CONDITION                                           
!===============================================================================================!
  dt     = abs(dtMax)
  time   = 0.0_mytype
  istart = 0
  istep  = 0
  CFL_new = CFL 

! initialize initial solution or read restart file
  if (readRestartFile.lt.0)  then
#if defined(BL)
  ! Blasius
  call initField_BL(part,x,z,rho,u,v,w,ien,pre,tem,mu,ka)
  rho_bl=rho;u_bl=u;v_bl=v;w_bl=w;ien_bl=ien;pre_bl=pre;tem_bl=tem;mu_bl=mu;ka_bl=ka
  !$acc update device(rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl)
#elif defined(CHA) 
  ! channel
  call initField_CHA(part,rho,u,v,w,ien,pre,tem,mu,ka)
#elif defined(TGV)
  ! TGV
  call initField_TGV(part,rho,u,v,w,ien,pre,tem,mu,ka)
#elif defined(DHC)
  ! TGV
  call initField_DHC(part,rho,u,v,w,ien,pre,tem,mu,ka)
#elif defined (ONED) 
  ! 1D wave
  call initField_1D(part,rho,u,v,w,ien,pre,tem,mu,ka)
#endif
! save restart at time t=0
    !$acc update host(rho,u,v,w,ien)
    call saveRestart(istart,time,dpdz,rho,u,v,w,ien,nHalo,part)
  else
! read restart file    
    istart = readRestartFile
    call loadRestart(istart,time,dpdz,rho,u,v,w,ien,nHalo,part)
    !$acc update device(rho,u,v,w,ien)
  endif

! update the other secondary variables
  if (nrank==0) write(stdout,*) 'updating initial pre,tem,mu,ka'
  call calcState_RE(rho,ien,pre,tem,mu,ka,1,xsize(1),1,xsize(2),1,xsize(3)) 
  call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
  if (nrank==0) write(stdout,*)
  if ((perBC(3) .eqv. .false.) .and. (BC_inl_rescale .eqv. .true.)) then
    ! rescaling inlet boundary
      call setBC_Inl_rescale(rho,u,v,w,ien,pre,tem,mu,ka)
  endif
  call print_pertBC()
! calculate vorticity
  call calcVort(vortx,vorty,vortz,strxz,u,v,w)

! write multiple planes  
  !$acc update host(rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz,strxz)
  if (yi_plane(1).gt.0)  then
    do i=1,size(yi_plane)
      call output2dPlane(part,nHalo,istart,2,yi_plane(i),rho,'ypl.','r.', u,'ypl.','u.', v,'ypl.','v.', w,'ypl.','w.', &
                         pre,'ypl.','p.', tem,'ypl.','t.', ien,'ypl.','e.', mu,'ypl.','mu.', &
                         ka,'ypl.','ka.',vortx,'ypl.','vortx.',vorty,'ypl.','vorty.',vortz,'ypl.','vortz.')
      enddo
  endif
  if (xi_plane(1).gt.0)  then  
    do i=1,size(xi_plane)
      if ((xi_plane(i) .eq. 1) .and. (perBC(1) .eqv. .false.)) then !planes at lower wall
        call output2dPlane(part,nHalo,istart,1,xi_plane(i),rho,'xpl.','r.', tem,'xpl.','t.', mu,'xpl.','mu.', &
                                                                                             strxz,'xpl.','strxz.')
      else
        call output2dPlane(part,nHalo,istart,1,xi_plane(i),rho,'xpl.','r.', u,'xpl.','u.', v,'xpl.','v.', w,'xpl.','w.', &
                           pre,'xpl.','p.', tem,'xpl.','t.', ien,'xpl.','e.', mu,'xpl.','mu.', &
                           ka,'xpl.','ka.',vortx,'xpl.','vortx.',vorty,'xpl.','vorty.',vortz,'xpl.','vortz.')
      endif
    enddo
  endif
  if (zi_plane(1).gt.0)  then 
    do i=1,size(zi_plane) 
      call output2dPlane(part,nHalo,istart,3,zi_plane(i),rho,'zpl.','r.', u,'zpl.','u.', v,'zpl.','v.', w,'zpl.','w.', &
                         pre,'zpl.','p.', tem,'zpl.','t.', ien,'zpl.','e.', mu,'zpl.','mu.', &
                         ka,'zpl.','ka.',vortx,'zpl.','vortx.',vorty,'zpl.','vorty.',vortz,'zpl.','vortz.')
      enddo
  endif

!===============================================================================================!
!                                          TIME ITERATION
!===============================================================================================!
  !$acc update device(rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz)
  wt_start = MPI_WTIME()
  wt_tmp = MPI_WTIME()

  ! calculate bulk properties for real time monitoring
  call cmpbulk(istart,wt_tmp,time,dt,CFL_new,rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz,velwb,dpdz)

  ! calculate timestep
  if (CFL.gt.0) call calcTimeStep(dt,CFL_new,rho,u,v,w,ien,mu,ka)

  ! loop over timesteps
  do istep=istart+1, istart+nsteps
    wt_tmp = MPI_WTIME()
    ! for channel flow: pressure gradient correction for unity bulk velocity
    if (dpdz /= 0.0) then 
      call cmpbulkvel(rho,w,velwb)
      dpdz = 0.999*dpdz - 0.5*(velwb - 1.0)
    endif
    if ((perBC(3) .eqv. .false.) .and. (BC_inl_rescale .eqv. .true.)) then
    ! rescaling inlet boundary
      call setBC_Inl_rescale(rho,u,v,w,ien,pre,tem,mu,ka)
    endif

    ! recalculate timestep if condition is met
    if ((CFL.gt.0).and.(mod(istep, intvCalcCFL).eq.0))  call calcTimeStep(dt,CFL_new,rho,u,v,w,ien,mu,ka)
    ! time integration

    call rk3(part,dt,istep,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time) !
    ! advance time
    time = time + dt

    ! calculate bulk properties for real time monitoring
    if (mod(istep,intvPrint).eq.0) then 
      call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
      call calcVort(vortx,vorty,vortz,strxz,u,v,w) 
      call cmpbulk(istep,wt_tmp,time,dt,CFL_new,rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz,velwb,dpdz)
    endif

    ! I/O planes if condition is met 
    if ((mod(istep-istart,intvSavePlanes).eq.0).and.(istep.ge.savePlanesAfter).and.(savePlanesAfter .ge. 0)) then 
      ! set boundary conditions
      call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
      call calcVort(vortx,vorty,vortz,strxz,u,v,w) 
      !$acc update host(rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz,strxz)
      ! write y-planes
      if (yi_plane(1).gt.0)  then
        do i=1,size(yi_plane)
          call output2dPlane(part,nHalo,istep,2,yi_plane(i),rho,'ypl.','r.', u,'ypl.','u.', v,'ypl.','v.', w,'ypl.','w.', &
                             pre,'ypl.','p.', tem,'ypl.','t.', ien,'ypl.','e.', mu,'ypl.','mu.', &
                             ka,'ypl.','ka.', vortx,'ypl.','vortx.',vorty,'ypl.','vorty.',vortz,'ypl.','vortz.')
        enddo
      endif
      ! write x-planes
      if (xi_plane(1).gt.0)  then
        do i=1,size(xi_plane)
          if ((xi_plane(i) .eq. 1) .and. (perBC(1) .eqv. .false.)) then !planes at lower wall
            call output2dPlane(part,nHalo,istep,1,xi_plane(i),rho,'xpl.','r.', tem,'xpl.','t.', mu,'xpl.','mu.', &
                                                                                                strxz,'xpl.','strxz.')
          else
            call output2dPlane(part,nHalo,istep,1,xi_plane(i),rho,'xpl.','r.', u,'xpl.','u.', v,'xpl.','v.', w,'xpl.','w.', &
                              pre,'xpl.','p.', tem,'xpl.','t.', ien,'xpl.','e.', mu,'xpl.','mu.', &
                              ka,'xpl.','ka.', vortx,'xpl.','vortx.',vorty,'xpl.','vorty.',vortz,'xpl.','vortz.')
          endif
        enddo
      endif
      ! write z-planes
      if (zi_plane(1).gt.0)  then
        do i=1,size(zi_plane)
          call output2dPlane(part,nHalo,istep,3,zi_plane(i),rho,'zpl.','r.', u,'zpl.','u.', v,'zpl.','v.', w,'zpl.','w.', &
                             pre,'zpl.','p.', tem,'zpl.','t.', ien,'zpl.','e.', mu,'zpl.','mu.', &
                             ka,'zpl.','ka.', vortx,'zpl.','vortx.',vorty,'zpl.','vorty.',vortz,'zpl.','vortz.')
        enddo
      endif
    endif

    ! I/O restart files if condition is met 
    if ((mod(istep-istart,intvSaveRestart) .eq. 0).and.(istep .ge. saveRestartAfter).and.(saveRestartAfter .ge. 0)) then
      call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
      !$acc update host(rho,u,v,w,ien)
      call saveRestart(istep,time,dpdz,rho,u,v,w,ien,nHalo,part)
    endif

    ! I/O statistics files if condition is met 
    if ((mod(istep-istart,intvSaveStats) .eq. 0).and.(istep .ge. saveStatsAfter).and.(saveStatsAfter .ge. 0)) then
      if (istep .eq. saveStatsAfter) then
        if (nrank==0) write(stdout,'(A, I10)') 'Begin of averaging at initial step:   ', istep
      else if ((istep .ge. saveStatsAfter) .and. (istep .ne. (istart+nsteps))) then
        if (nrank==0) write(stdout,'(A, I10)') 'next averaging at step:               ', istep
      else if ((istep .ge. saveStatsAfter) .and. (istep .eq. (istart+nsteps))) then
        if (nrank==0) write(stdout,'(A, I10)') 'final averaging at step:              ', istep
      endif
      call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
      call calcStats(qave,qtime,factAvg,countAvg,rho,u,v,w,ien,pre,tem,mu,ka)
      if (istep .eq. (istart+nsteps)) then
        !$acc update host(rho,u,v,w,ien,pre,tem,mu,ka,qave,qtime)
        call saveStats(part,istep,dt,qave,qtime,factAvg,countAvg)
      endif
    endif

    ! I/O file for real-time parameter modifications
    if ((mod(istep,intvReadParam) .eq. 0) .and. (intvReadParam .gt. 0)) call io_readParams_variation()
    ! simulation can be stopped with *abortSimulation*
    if (abortSimulation .eqv. .true.) then 
      !$acc update host(rho,u,v,w,ien)
      call saveRestart(istep,time,dpdz,rho,u,v,w,ien,nHalo,part)    
      call decomp_2d_finalize
      call mpi_finalize(ierr)
      stop
    endif
  enddo
  if (nrank==0) then
    write(*,*)
    write (stdout, *) "o--------------------------------------------------o"
    write(*,*) 'SIMULATE done!'
  endif
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
  deallocate(vortx)
  deallocate(vorty)
  deallocate(vortz)
  deallocate(strxz)
  deallocate(qave)
  deallocate(qtime)
#if defined(BL)
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
#endif
  ! stop the mpi 
  call decomp_2d_finalize
  call mpi_finalize(ierr)


end program
