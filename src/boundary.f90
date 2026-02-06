! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! boundary conditions module

module mod_boundary
  use mod_param
  use mod_grid
  use mod_halo
  use mod_math
  use mod_timer
  implicit none
  real(mytype), allocatable, dimension(:,:) :: dudt_pert
  real(mytype), allocatable, dimension(:,:) :: w_LST_real, u_LST_real, tem_LST_real
  real(mytype), allocatable, dimension(:,:) :: w_LST_imag, u_LST_imag, tem_LST_imag
  real(mytype), allocatable, dimension(:) :: au_rcy,av_rcy,aw_rcy,ap_rcy,at_rcy
  real(mytype), allocatable, dimension(:) :: x_rcy_inner,x_rcy_outer,weight_func 
  real(mytype) :: beta_rescale,delta_rcy,A_FF
  integer :: nLST
contains


  ! initialization boundary condition module
  subroutine init_BC()
    use decomp_2d
    use mod_param
    use mod_eos
    use mod_halo
    use mod_finitediff
    implicit none
    integer                            :: i, ierr, fs                     ! Loop index, ierr
    logical                            :: bcexists                        ! Existence check for (physical) boundary conditions
    ! prescribing eigenfunctions based on Linear Stability Theory
    if (BC_inl == "inlet_lst") then
      ! loading of the eigenfunctions mesh
      inquire(FILE='initBL/inputDNS/prof_x_LST.bin', SIZE=fs); nLST = fs/8
      if (fs < 0) then
        write(stdout,*) "ERROR! no input LST file for inlet_lst"
        call decomp_2d_finalize
        call mpi_finalize(ierr) 
        stop
      endif
      allocate(w_LST_real(1:xsize(1),1:xsize(2)),u_LST_real(1:xsize(1),1:xsize(2)),tem_LST_real(1:xsize(1),1:xsize(2)))
      allocate(w_LST_imag(1:xsize(1),1:xsize(2)),u_LST_imag(1:xsize(1),1:xsize(2)),tem_LST_imag(1:xsize(1),1:xsize(2)))
      !$acc enter data create(w_LST_real,u_LST_real,tem_LST_real,w_LST_imag,u_LST_imag,tem_LST_imag)
      ! initialization of the eigenfunctions
      call initLST(nLST,w_LST_real,u_LST_real,tem_LST_real,w_LST_imag,u_LST_imag,tem_LST_imag)
    endif
    ! recycling and rescaling
    if (BC_inl_rescale .eqv. .true.)then
      call init_rescale()
    endif
    if (perBC(1) .eqv. .true.) then
      BC_top = 'periodic'
      BC_bot = 'periodic'
    endif
    if (perBC(3) .eqv. .true.) then
      BC_inl = 'periodic'
      BC_out = 'periodic'
    endif
    if (perBC(2) .eqv. .true.) then
      BC_span = 'periodic'
    endif

    ! print the boundary conditions
    if (nrank == 0) then
      write(stdout,* ) 'Boundary conditions'
      write(stdout,'(A)') 'o--------------------------------------------------o'
      write(stdout,'(A, F10.4)') 'Boundary conditions                        done!'
#if defined(BL) 
      if (BC_bot(1:5)==wall_bc) then
        write(stdout,* ) 'Correct BC_bot in initBL and config.h'    
      else
        write(stdout,* ) 'Mismatch BC_bot between initBL and config.h, check both again!'
        call decomp_2d_finalize
        call mpi_finalize(ierr) 
        stop
      endif
#endif
    write(stdout,'(A, A10)') 'Top boundary:                         ',BC_top
    write(stdout,'(A, A10)') 'Bottom boundary:                      ',BC_bot
#if defined(BL)
      write(stdout,'(A, F10.4)') 'T_wall (T/T_inf):                     ',Twall_bot
#elif defined(CHA)
      write(stdout,'(A, F10.4)') 'Top T_wall (T/T_inf):                 ',Twall_top
      write(stdout,'(A, F10.4)') 'Bottom T_wall (T/T_inf):              ',Twall_bot
#elif defined(DHC)
      write(stdout,'(A, F10.4)') 'Inlet T_wall (T/T_ref):               ',Twall_inl
      write(stdout,'(A, F10.4)') 'Outlet T_wall (T/T_ref):              ',Twall_out
#endif
    write(stdout,'(A, A10)') 'Inlet boundary:                       ',BC_inl
    write(stdout,'(A, A11)') 'Outlet boundary:                      ',BC_out
    write(stdout,'(A, A10)') 'Spanwise boundary:                      ',BC_span
    write(stdout,'(A)') 'o--------------------------------------------------o'
    write(stdout,* )
    endif
    if ((p_row == 1) .and. (nrank == i)) then
      if ((neigh%inlet)) then 
         write(stdout,'(A, I10)') 'Inlet boundary, proc=                ',nrank
      endif
    endif
    if ((p_row == 1) .and. (nrank==p_col-1)) then
      if ((neigh%outlet))  then 
         write(stdout,'(A, I10)')'Outlet boundary, proc=                ',nrank
         write(stdout,* )
      endif
    endif
  end subroutine


! initialization of the eigenfunctions
  subroutine initLST(nLST,w_LST_real,u_LST_real,tem_LST_real,w_LST_imag,u_LST_imag,tem_LST_imag)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_eos_var
  use mod_grid
  use mod_halo
  use mod_math
  implicit none
  integer :: i,j,k, fs, ierr, kk, lenr, kBC, nLST
  complex(mytype), allocatable, dimension(:) :: wRead_fluc, uRead_fluc, tRead_fluc
  real(mytype), allocatable, dimension(:) :: wRead, uRead, tRead
  real(mytype), allocatable, dimension(:) :: w_fluct_real, u_fluct_real, tem_fluct_real, w_fluct_imag, u_fluct_imag, tem_fluct_imag
  real(mytype), allocatable, dimension(:) :: wIntp_real, uIntp_real, tIntp_real, wIntp_imag, uIntp_imag, tIntp_imag
  real(mytype), allocatable, dimension(:) :: xLST
  real(mytype) :: dummy, scaling_delta, x0, fact, time !!!
  real(mytype), dimension(:,:) :: w_LST_real, u_LST_real, tem_LST_real, w_LST_imag, u_LST_imag, tem_LST_imag
  if (nrank==0) then
    write(stdout,*) 'Initializing LST eigenfunction'
  endif
  fs = nLST * 8
  ! allocate variables
  allocate(xLST(nLST))      
  allocate(wRead(nLST*2)); allocate(wRead_fluc(nLST)) 
  allocate(uRead(nLST*2)); allocate(uRead_fluc(nLST))
  allocate(tRead(nLST*2)); allocate(tRead_fluc(nLST))
  allocate(w_fluct_real(nLST),u_fluct_real(nLST),tem_fluct_real(nLST))
  allocate(w_fluct_imag(nLST),u_fluct_imag(nLST),tem_fluct_imag(nLST))
  allocate(wIntp_real(nLST),uIntp_real(nLST),tIntp_real(nLST))
  allocate(wIntp_imag(nLST),uIntp_imag(nLST),tIntp_imag(nLST))
  ! open eigenfunction files
  open(9,file='initBL/inputDNS/prof_x_LST.bin',access='direct',recl=fs); read(9,rec=1) (xLST(i),i=1,nLST);close(9)
  open(9,file='initBL/inputDNS/prof_w_LST.bin',access='direct',recl=fs*2); read(9,rec=1) (wRead(i),i=1,nLST*2);close(9)
  open(9,file='initBL/inputDNS/prof_u_LST.bin',access='direct',recl=fs*2); read(9,rec=1) (uRead(i),i=1,nLST*2);close(9)
  open(9,file='initBL/inputDNS/prof_T_LST.bin',access='direct',recl=fs*2); read(9,rec=1) (tRead(i),i=1,nLST*2);close(9)
  do i=1,nLST
    wRead_fluc(i)=cmplx(wRead(i),wRead(nLST+i))
    uRead_fluc(i)=cmplx(uRead(i),uRead(nLST+i))
    tRead_fluc(i)=cmplx(tRead(i),tRead(nLST+i))
  enddo
  do i=1,nLST
    w_fluct_real(i)=real(wRead_fluc(i))
    u_fluct_real(i)=real(uRead_fluc(i))
    tem_fluct_real(i)=real(tRead_fluc(i))

    w_fluct_imag(i)=aimag(wRead_fluc(i))
    u_fluct_imag(i)=aimag(uRead_fluc(i))
    tem_fluct_imag(i)=aimag(tRead_fluc(i))
  enddo
  ! interpolation
  call spline(xLST, w_fluct_real, nLST, wIntp_real)
  call spline(xLST, u_fluct_real, nLST, uIntp_real)
  call spline(xLST, tem_fluct_real, nLST, tIntp_real)
  call spline(xLST, w_fluct_imag, nLST, wIntp_imag)
  call spline(xLST, u_fluct_imag, nLST, uIntp_imag)
  call spline(xLST, tem_fluct_imag, nLST, tIntp_imag)
    do j=1,xsize(2)
      do i=1,xsize(1)
        call splint(xLST, w_fluct_real, wIntp_real, nLST, x(i),  w_LST_real(i,j))
        call splint(xLST, u_fluct_real, uIntp_real, nLST, x(i),  u_LST_real(i,j))
        call splint(xLST, tem_fluct_real, tIntp_real, nLST, x(i),  tem_LST_real(i,j))
        call splint(xLST, w_fluct_imag, wIntp_imag, nLST, x(i),  w_LST_imag(i,j))
        call splint(xLST, u_fluct_imag, uIntp_imag, nLST, x(i),  u_LST_imag(i,j))
        call splint(xLST, tem_fluct_imag, tIntp_imag, nLST, x(i),  tem_LST_imag(i,j))
      enddo
    enddo
  !$acc update device(w_LST_real,u_LST_real,tem_LST_real,w_LST_imag,u_LST_imag,tem_LST_imag)
  ! deallocate variables
  deallocate(wRead); deallocate(wRead_fluc)
  deallocate(uRead); deallocate(uRead_fluc)
  deallocate(tRead); deallocate(tRead_fluc)
  deallocate(w_fluct_real,u_fluct_real,tem_fluct_real)
  deallocate(w_fluct_imag,u_fluct_imag,tem_fluct_imag)
  deallocate(wIntp_real,uIntp_real,tIntp_real)
  deallocate(wIntp_imag,uIntp_imag,tIntp_imag)
end subroutine


! initialization of the rescale for turbulent boundary layer
subroutine init_rescale()
  use decomp_2d
  use mod_param
  use mod_grid
  use mod_halo
  implicit none
  real(mytype) :: xdelta_inl
  integer :: i,ierr
  real(mytype), parameter :: prt   = 0.89_mytype ! Turbulent Prandtl number

  allocate(au_rcy(1:xsize(1)),av_rcy(1:xsize(1)),aw_rcy(1:xsize(1)),ap_rcy(1:xsize(1)),at_rcy(1:xsize(1)))
  allocate(x_rcy_inner(1:xsize(1)),x_rcy_outer(1:xsize(1)),weight_func(1:xsize(1)))
  ! reading mean profiles from files
  open(11,file = 'preproc/turbRR/mean_values.txt',form='formatted',status="old",action="read")
  read(11,*)
  do i=1,xsize(1)
    read(11,*) au_rcy(i),av_rcy(i),aw_rcy(i),ap_rcy(i),at_rcy(i)
  enddo
  close(11)
  ! obtaining rescaling parameter 'beta'
  do i=1,xsize(1)
    ! 99% free-stream velocity
    if (aw_rcy(i) .gt. 0.99_mytype) then 
      delta_rcy = x(i)
      exit
    endif
  enddo
  if (delta_rcy .le. delta_inl) then
    write(stdout,* ) 'ERROR: delta_inl is too large. please set the smaller value than:', delta_rcy
    call decomp_2d_finalize
    call mpi_finalize(ierr)
    stop
  endif
  beta_rescale = (delta_rcy/delta_inl)**0.1
  ! setting weighting function (ideal gas)
  A_FF = dsqrt(((ig_gam-1.0_mytype)*0.5_mytype*(Ma)**2*prt)/ &
       (1.0_mytype + (ig_gam-1.0_mytype)*0.5_mytype*(Ma)**2*prt))
  x_rcy_inner = x*beta_rescale
  x_rcy_outer = x*beta_rescale**10.0_mytype
  do i=1,xsize(1)
    xdelta_inl = x(i)/delta_inl
    weight_func(i) = 0.5_mytype*(1.0_mytype + (tanh(4.0_mytype*(xdelta_inl-0.2_mytype)/ &
                   ((1.0_mytype - 2.0_mytype*0.2_mytype)*xdelta_inl+0.2_mytype))/tanh(4.0_mytype)))
    if (xdelta_inl .ge. 1.0_mytype) then
      weight_func(i) = 1.0_mytype
    endif
  enddo
  if (nrank == 0) then
    write(stdout,* ) 'Recycling and rescaling parameters'
    write(stdout,'(A)') '------------------------------------------------'
    write(stdout,'(A, F10.4)') 'Recycling position:                   ',z_recycle
    write(stdout,'(A, F10.4)') 'Recycling Reynolds number:            ',(Re)**0.5_mytype*(zStartDNS+z_recycle)**0.5_mytype
    write(stdout,'(A, F10.4)') 'Beta parameter:                       ',beta_rescale
    write(stdout,* )
  endif
  !$acc enter data copyin(au_rcy,av_rcy,aw_rcy,ap_rcy,at_rcy) 
  !$acc enter data copyin(x_rcy_inner,x_rcy_outer,xdelta_inl,weight_func) 
end subroutine


! set all boundary conditions: bottom, top, inlet, outlet
subroutine setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_halo
  use mod_finitediff
  implicit none
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl
  real(mytype) :: time
  TYPE (DECOMP_INFO), intent(IN) :: part

  ! communication: halo cells data transfer
#if defined(_OPENACC)
  ! GPU
  call haloUpdateMult_i(perBC,xsize,rho,u,v,w,ien,pre,tem,mu,ka)  
  call haloUpdate9_jk(perBC,xsize,rho,u,v,w,ien,pre,tem,mu,ka) 
#else
  ! CPU
  call haloUpdateMult_CPU(perBC,xsize,rho,u,v,w,ien,pre,tem,mu,ka)  
#endif
  ! set halo cells for boundary conditions
  ! in case perBC is set to .true., periodic boundary conditions are applied
  if (perBC(1) .eqv. .false.) then 
    ! bottom boundary
    call setBC_Bot(rho,u,v,w,ien,pre,tem,mu,ka,time)
    ! top boundary
    call setBC_Top(rho,u,v,w,ien,pre,tem,mu,ka)
  endif
  if ((perBC(3) .eqv. .false.) .and. (BC_inl_rescale .eqv. .false.) .and. (neigh%inlet)) then 
    ! standard inlet boundary 
    call setBC_Inl(rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
  endif
  if ((perBC(3) .eqv. .false.) .and. (neigh%outlet)) then 
    ! outlet boundary
    call setBC_Out(rho,u,v,w,ien,pre,tem,mu,ka)
  endif
end subroutine


! halo cells: bottom boundary
subroutine setBC_Bot(rho,u,v,w,ien,pre,tem,mu,ka,time)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_halo
  use mod_finitediff
  use mod_perturbation
  implicit none
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(mytype) :: time
  integer :: j,k,c, ierr, jm,km, iBC
  ! at the first mesh cell in x-direction
  iBC = 1
  jm = xsize(2)
  km = xsize(3)
  !$acc parallel loop collapse(2) default(present) async(1) 
  do k=1,km
    do j=1,jm
      u(iBC,j,k) = 0.0_mytype
      v(iBC,j,k) = 0.0_mytype
      w(iBC,j,k) = 0.0_mytype
      dudt_pert(j,k) = 0.0_mytype 
    enddo
  enddo
  !$acc wait(1)
  ! adiabatic standard boundary condition
  if (BC_bot == "adiab_std") then 
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm
        pre(iBC,j,k) = 0.0_mytype
        tem(iBC,j,k) = 0.0_mytype
      enddo
    enddo
    !$acc wait(1)    
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm
        !$acc loop seq
        do c = 1,2
          pre(iBC,j,k) = pre(iBC,j,k) - c1_FD2(c)*pre(iBC+c,j,k)/c1_FD2(0)
          tem(iBC,j,k) = tem(iBC,j,k) - c1_FD2(c)*tem(iBC+c,j,k)/c1_FD2(0)
        enddo
      enddo
    enddo
    !$acc wait(1)
    call calcState_PT(pre,tem,rho,ien,mu,ka,iBC,iBC,1,jm,1,km)
    ! update halo cells
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm 
        do c=1,nHalo 
          rho( iBC-c,j,k) = 2.0_mytype*rho(iBC,j,k) - rho(iBC+c,j,k)
            u( iBC-c,j,k) = 2.0_mytype*  u(iBC,j,k) -   u(iBC+c,j,k)
            v( iBC-c,j,k) = 2.0_mytype*  v(iBC,j,k) -   v(iBC+c,j,k)
            w( iBC-c,j,k) = 2.0_mytype*  w(iBC,j,k) -   w(iBC+c,j,k)
          ien( iBC-c,j,k) = 2.0_mytype*ien(iBC,j,k) - ien(iBC+c,j,k)
          pre( iBC-c,j,k) =                           pre(iBC+c,j,k)
          tem( iBC-c,j,k) =                           tem(iBC+c,j,k)
           mu( iBC-c,j,k) = 2.0_mytype* mu(iBC,j,k) -  mu(iBC+c,j,k)
           ka( iBC-c,j,k) = 2.0_mytype* ka(iBC,j,k) -  ka(iBC+c,j,k)
        enddo
      enddo
    enddo
  ! adiabatic non-reflecting boundary condition
  else if (BC_bot == "adiab_nrbc") then 
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm 
        do c=1,nHalo 
          rho( iBC-c,j,k) = 2.0_mytype*rho(iBC,j,k) - rho(iBC+c,j,k)
            u( iBC-c,j,k) = 2.0_mytype*  u(iBC,j,k) -   u(iBC+c,j,k)
            v( iBC-c,j,k) = 2.0_mytype*  v(iBC,j,k) -   v(iBC+c,j,k)
            w( iBC-c,j,k) = 2.0_mytype*  w(iBC,j,k) -   w(iBC+c,j,k)
          ien( iBC-c,j,k) = 2.0_mytype*ien(iBC,j,k) - ien(iBC+c,j,k)
          pre( iBC-c,j,k) =                           pre(iBC+c,j,k)
          tem( iBC-c,j,k) =                           tem(iBC+c,j,k)
           mu( iBC-c,j,k) = 2.0_mytype* mu(iBC,j,k) -  mu(iBC+c,j,k)
           ka( iBC-c,j,k) = 2.0_mytype* ka(iBC,j,k) -  ka(iBC+c,j,k)
        enddo
      enddo
    enddo
  ! isothermal standard boundary condition
  else if ( BC_bot == "isoth_std" ) then 
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm
        pre(iBC,j,k) = 0.0_mytype
        ! prescribing the wall temperature
        tem(iBC,j,k) = Twall_bot
      enddo
    enddo
    !$acc wait(1)    
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm
        !$acc loop seq
        do c = 1,2
          pre(iBC,j,k) = pre(iBC,j,k) - c1_FD2(c)*pre(iBC+c,j,k)/c1_FD2(0)
        enddo
      enddo
    enddo
    !$acc wait(1)    
    call calcState_PT(pre,tem,rho,ien,mu,ka, iBC,iBC,1,jm,1,km)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm 
        do c=1,nHalo 
          rho( iBC-c,j,k) = 2.0_mytype*rho(iBC,j,k) - rho(iBC+c,j,k)
            u( iBC-c,j,k) = 2.0_mytype*  u(iBC,j,k) -   u(iBC+c,j,k)
            v( iBC-c,j,k) = 2.0_mytype*  v(iBC,j,k) -   v(iBC+c,j,k)
            w( iBC-c,j,k) = 2.0_mytype*  w(iBC,j,k) -   w(iBC+c,j,k)
          ien( iBC-c,j,k) = 2.0_mytype*ien(iBC,j,k) - ien(iBC+c,j,k)
          pre( iBC-c,j,k) =                           pre(iBC+c,j,k)
          tem( iBC-c,j,k) = 2.0_mytype*tem(iBC,j,k) - tem(iBC+c,j,k)
           mu( iBC-c,j,k) = 2.0_mytype* mu(iBC,j,k) -  mu(iBC+c,j,k)
           ka( iBC-c,j,k) = 2.0_mytype* ka(iBC,j,k) -  ka(iBC+c,j,k)
        enddo
      enddo
    enddo
 ! isothermal non-reflecting boundary condition
  else if (BC_bot == "isoth_nrbc") then 
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm 
        ! prescribing the wall temperature
        tem(iBC,j,k) = Twall_bot
      enddo
    enddo
    !$acc wait(1)
    call calcState_rT(rho,tem,ien,pre,mu,ka, iBC,iBC,1,jm,1,km)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm 
        do c=1,nHalo 
          rho( iBC-c,j,k) = 2.0_mytype*rho(iBC,j,k) - rho(iBC+c,j,k)
            u( iBC-c,j,k) = 2.0_mytype*  u(iBC,j,k) -   u(iBC+c,j,k)
            v( iBC-c,j,k) = 2.0_mytype*  v(iBC,j,k) -   v(iBC+c,j,k)
            w( iBC-c,j,k) = 2.0_mytype*  w(iBC,j,k) -   w(iBC+c,j,k)
          ien( iBC-c,j,k) = 2.0_mytype*ien(iBC,j,k) - ien(iBC+c,j,k)
          pre( iBC-c,j,k) =                           pre(iBC+c,j,k)
          tem( iBC-c,j,k) = 2.0_mytype*tem(iBC,j,k) - tem(iBC+c,j,k)
           mu( iBC-c,j,k) = 2.0_mytype* mu(iBC,j,k) -  mu(iBC+c,j,k)
           ka( iBC-c,j,k) = 2.0_mytype* ka(iBC,j,k) -  ka(iBC+c,j,k)
        enddo
      enddo
    enddo
  ! non-existing boundary condition  
  else
    if (nrank .eq. 0) write(stdout,*) "unknown BC_bottom: ", BC_bot
    call decomp_2d_finalize
    call mpi_finalize(ierr)
    stop
  endif
  !$acc wait(1)
  ! prescribing the disturbance strip
  if (pert_calc == 1) then
    ! comment and uncomment the desired perturbation
    !call perturbationFrankoLele(u,time, 0,dudt_pert)
    call perturbationSayadi(u,time,0,dudt_pert)
  ! prescribing eigenfunctions at the inlet 
  elseif ((pert_calc == 1) .AND. (BC_inl == "inlet_lst")) then
    if  (nrank .eq. 0) then
      write(stdout,*) "inlet_lst:", BC_inl
      write(stdout,*) "ERROR! with inlet_lst pert_calc should be set to 0"
    endif
    call decomp_2d_finalize
    call mpi_finalize(ierr) 
    stop
  endif
end subroutine


! RHS: bottom boundary 
subroutine setBC_RHS_Bot(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,rho,u,v,w,ien,pre)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_halo
  use mod_perturbation
  use mod_finitediff
  implicit none
  real(mytype), dimension(:,:,:) :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre
  real(mytype), dimension(5) :: d,L
  real(mytype) :: Kfact, sigm, sos, fac, dp, drho, du, dv, dw, ht
  real(mytype) :: xa, rhoa, ua, iena
  integer :: iBC,j,k,jm,km,c, ierr, k1
  ! at the first mesh cell in x-direction
  iBC = 1
  jm = xsize(2)
  km = xsize(3)
  ! calculate characteristics for adiabatic non-reflecting boundary condition
  if (BC_bot == "adiab_nrbc") then 
    ! if perturbation is active
    if (pert_calc == 1) then
        !$acc parallel loop collapse(2) default(present) private(sos,fac,dp,drho,du,d,L,ht) async(1)
         do k = 1,km
            do j = 1,jm
              xa = xp(iBC)
              rhoa = rho(iBC,j,k)
              ua = u(iBC,j,k)
              iena = ien(iBC,j,k)
              call calcSOS(rhoa,iena,sos)
              call calcFac(rhoa,iena,fac)
              call calc_conv_FD_ddx(dp,pre,iBC,j,k,nHalo,xa,dx)
              call calc_conv_FD_ddx(drho,rho,iBC,j,k,nHalo,xa,dx)
              call calc_conv_FD_ddx(du,u,iBC,j,k,nHalo,xa,dx)
              ! amplitudes of characteristic waves
              L(1) = (ua - sos)*(dp - rhoa*sos*du)
              L(2) = ua*(dp - (sos**2)*drho)
              L(3) = 0.0_mytype
              L(4) = 0.0_mytype
              L(5) = L(1) - 2*sos*(rhoa*dudt_pert(j,k) + Ri_unit*(rhoa-1.0_mytype))
              ! time variation of the wave amplitudes
              d(1) = (0.5_mytype*(L(5) + L(1)) - L(2))/sos**2
              d(2) = (L(5) - L(1))/(2*rhoa*sos)
              d(3) = L(3)
              d(4) = L(4)
              d(5) = fac*L(2)/sos**2
              ht = iena + pre(iBC,j,k)/rhoa + 0.5_mytype*(ua**2 + v(iBC,j,k)**2 + w(iBC,j,k)**2)
              ! correct the right hand side
              rhs_r(iBC,j,k) = rhs_r(iBC,j,k) - d(1)
              rhs_u(iBC,j,k) = ua * rhs_r(iBC,j,k) + dudt_pert(j,k) * rhoa &
                               - d(1)* ua - d(2)*rhoa
              rhs_v(iBC,j,k) = 0.0_mytype
              rhs_w(iBC,j,k) = 0.0_mytype
              rhs_e(iBC,j,k) = rhs_e(iBC,j,k) - d(1)*ht &
                                              - d(2)*rhoa*ua & 
                                              - d(5)                                            
            enddo
          enddo
    ! if perturbation is not active
    else
      !$acc parallel loop collapse(2) default(present) private(sos,fac,dp,du,d,L,ht) async(1)
      do k = 1,km
        do j = 1,jm
        xa = xp(iBC)
        ua = u(iBC,j,k)
        rhoa = rho(iBC,j,k)
        iena = ien(iBC,j,k)
        call calcSOS(rhoa,iena,sos)
        call calcFac(rhoa,iena,fac)
        call calc_conv_FD_ddx(dp,pre,iBC,j,k,nHalo,xa,dx)
        call calc_conv_FD_ddx(du,u,iBC,j,k,nHalo,xa,dx)
        ! amplitudes of characteristic waves
        L(1) = (ua - sos)*(dp - rhoa*sos*du)
        L(2) = 0.0_mytype
        L(3) = 0.0_mytype
        L(4) = 0.0_mytype
        L(5) = L(1) - 2.0*sos*Ri_unit*(rhoa-1.0_mytype)
        ! time variation of the wave amplitudes
        d(1) = (0.5_mytype*(L(5) + L(1)) - L(2))/sos**2
        d(2) = (L(5) - L(1))/(2*rhoa*sos)
        d(3) = L(3)
        d(4) = L(4)
        d(5) = fac*L(2)/sos**2
        ht = iena + pre(iBC,j,k)/rhoa + 0.5_mytype*(ua**2 + v(iBC,j,k)**2 + w(iBC,j,k)**2)
        ! correct the right hand side
        rhs_r(iBC,j,k) = rhs_r(iBC,j,k) - d(1)
        rhs_u(iBC,j,k) = 0.0_mytype
        rhs_v(iBC,j,k) = 0.0_mytype
        rhs_w(iBC,j,k) = 0.0_mytype
        rhs_e(iBC,j,k) = rhs_e(iBC,j,k) - d(1)*ht &
                                        - d(5)
        enddo
      enddo
    endif
  ! calculate characteristics for adiabatic non-reflecting boundary condition
  else if (BC_bot == "isoth_nrbc") then 
    ! if perturbation is active
    if (pert_calc==1) then
      !$acc parallel loop collapse(2) default(present) private(sos,fac,dp,drho,du,d,L,ht) async(1)
      do k = 1,km
        do j = 1,jm
          xa = xp(iBC)
          rhoa = rho(iBC,j,k)
          ua = u(iBC,j,k)
          iena = ien(iBC,j,k)
          call calcSOS(rhoa,iena,sos)
          call calcFac(rhoa,iena,fac)
          call calc_conv_FD_ddx(dp,pre,iBC,j,k,nHalo,xa,dx)
          call calc_conv_FD_ddx(drho,rho,iBC,j,k,nHalo,xa,dx)
          call calc_conv_FD_ddx(du,u,iBC,j,k,nHalo,xa,dx)
          ! amplitudes of characteristic waves
          L(1) = (ua - sos)*(dp - rhoa*sos*du)
          L(2) = ua*(dp - (sos**2)*drho)
          L(3) = 0.0_mytype
          L(4) = 0.0_mytype 
          L(5) = L(1) - 2*sos*(rhoa*dudt_pert(j,k) + Ri_unit*(rhoa-1.0_mytype))
          ! time variation of the wave amplitudes
          d(1) = (0.5_mytype*(L(5) + L(1)) - L(2))/sos**2
          d(2) = (L(5) - L(1))/(2*rhoa*sos)
          d(3) = L(3)
          d(4) = L(4)
          d(5) = fac*L(2)/sos**2
          ! correct the right hand side
          rhs_r(iBC,j,k) = rhs_r(iBC,j,k) - d(1)
          rhs_u(iBC,j,k) = ua * rhs_r(iBC,j,k) + dudt_pert(j,k) * rhoa &
                           - d(1)* ua - d(2)*rhoa
          rhs_v(iBC,j,k) = 0.0_mytype
          rhs_w(iBC,j,k) = 0.0_mytype
          rhs_e(iBC,j,k) = 0.0_mytype
        enddo
      enddo
    ! if perturbation is not active
    else
      !$acc parallel loop collapse(2) default(present) private(sos,fac,dp,du,d,L,ht) async(1)
      do k = 1,km
        do j = 1,jm
          xa = xp(iBC)
          rhoa = rho(iBC,j,k)
          iena = ien(iBC,j,k)
          call calcSOS(rhoa,iena,sos)
          call calcFac(rhoa,iena,fac)
          call calc_conv_FD_ddx(dp,pre,iBC,j,k,nHalo,xa,dx)
          call calc_conv_FD_ddx(du,u,iBC,j,k,nHalo,xa,dx)
          ! amplitudes of characteristic waves
          L(1) = (u(iBC,j,k) - sos)*(dp - rhoa*sos*du)
          L(2) = 0.0_mytype
          L(3) = 0.0_mytype
          L(4) = 0.0_mytype
          L(5) = L(1) - 2*sos*Ri_unit*(rhoa-1.0_mytype)
          ! time variation of the wave amplitudes
          d(1) = (0.5_mytype*(L(5) + L(1)) - L(2))/sos**2
          d(2) = (L(5) - L(1))/(2*rhoa*sos)
          d(3) = L(3)
          d(4) = L(4)
          d(5) = fac*L(2)/sos**2
          ! correct the right hand side
          rhs_r(iBC,j,k) = rhs_r(iBC,j,k) - d(1)
          rhs_u(iBC,j,k) = 0.0_mytype
          rhs_v(iBC,j,k) = 0.0_mytype
          rhs_w(iBC,j,k) = 0.0_mytype
          rhs_e(iBC,j,k) = 0.0_mytype
        enddo
      enddo
    endif
  else
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k=1,km
      do j=1,jm
        rhs_r(iBC,j,k) = 0.0_mytype
        rhs_u(iBC,j,k) = 0.0_mytype
        rhs_v(iBC,j,k) = 0.0_mytype
        rhs_w(iBC,j,k) = 0.0_mytype
        rhs_e(iBC,j,k) = 0.0_mytype
      enddo
    enddo
  endif
end subroutine


! halo cells: top boundary
subroutine setBC_Top(rho,u,v,w,ien,pre,tem,mu,ka)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_halo
  use mod_finitediff
  implicit none
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(mytype), dimension(5) :: d,L
  real(mytype) :: Kfact, sigm, sos, fac, dp, drho, du, dv, dw
  integer :: iBC, ierr,c, jm, km, j, k 
  ! at the last mesh cell in x-direction
  iBC = xsize(1)
  jm = xsize(2)
  km = xsize(3)
  ! freestream boundary condition
  if (BC_top == "free_nrbc") then 
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm 
        do c=1,nHalo 
          rho(iBC+c,j,k) =  2.0_mytype*rho(iBC,j,k) - rho(iBC-c,j,k)
            u(iBC+c,j,k) =  2.0_mytype*  u(iBC,j,k) -   u(iBC-c,j,k)
            v(iBC+c,j,k) =  2.0_mytype*  v(iBC,j,k) -   v(iBC-c,j,k)
            w(iBC+c,j,k) =  2.0_mytype*  w(iBC,j,k) -   w(iBC-c,j,k)
          ien(iBC+c,j,k) =  2.0_mytype*ien(iBC,j,k) - ien(iBC-c,j,k)
          pre(iBC+c,j,k) =  2.0_mytype*pre(iBC,j,k) - pre(iBC-c,j,k)
          tem(iBC+c,j,k) =  2.0_mytype*tem(iBC,j,k) - tem(iBC-c,j,k)
           mu(iBC+c,j,k) =  2.0_mytype* mu(iBC,j,k) -  mu(iBC-c,j,k)
           ka(iBC+c,j,k) =  2.0_mytype* ka(iBC,j,k) -  ka(iBC-c,j,k)
        enddo
      enddo
    enddo
  ! adiabatic standard boundary condition
  else if (BC_top == "adiab_std") then
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k=1,km
      do j=1,jm
        u(iBC,j,k) = 0.0_mytype
        v(iBC,j,k) = 0.0_mytype
        w(iBC,j,k) = 0.0_mytype
      enddo
    enddo
    !$acc wait(1)
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm
        pre(iBC,j,k) = 0.0_mytype
        tem(iBC,j,k) = 0.0_mytype
      enddo
    enddo
    !$acc wait(1)    
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm
        !$acc loop seq
        do c = -2,-1
          pre(iBC,j,k) = pre(iBC,j,k) - c1_BD2(c)*pre(iBC+c,j,k)/c1_BD2(0)
          tem(iBC,j,k) = tem(iBC,j,k) - c1_BD2(c)*tem(iBC+c,j,k)/c1_BD2(0)
        enddo
      enddo
    enddo
    !$acc wait(1)
    call calcState_PT(pre,tem,rho,ien,mu,ka,iBC,iBC,1,jm,1,km)
    ! update halo cells
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm 
        do c=1,nHalo 
          rho( iBC+c,j,k) = 2.0_mytype*rho(iBC,j,k) - rho(iBC-c,j,k)
            u( iBC+c,j,k) = 2.0_mytype*  u(iBC,j,k) -   u(iBC-c,j,k)
            v( iBC+c,j,k) = 2.0_mytype*  v(iBC,j,k) -   v(iBC-c,j,k)
            w( iBC+c,j,k) = 2.0_mytype*  w(iBC,j,k) -   w(iBC-c,j,k)
          ien( iBC+c,j,k) = 2.0_mytype*ien(iBC,j,k) - ien(iBC-c,j,k)
          pre( iBC+c,j,k) =                           pre(iBC-c,j,k)
          tem( iBC+c,j,k) =                           tem(iBC-c,j,k)
           mu( iBC+c,j,k) = 2.0_mytype* mu(iBC,j,k) -  mu(iBC-c,j,k)
           ka( iBC+c,j,k) = 2.0_mytype* ka(iBC,j,k) -  ka(iBC-c,j,k)
        enddo
      enddo
    enddo
  ! adiabatic non-reflecting boundary condition
  else if (BC_top == "adiab_nrbc") then
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k=1,km
      do j=1,jm
        u(iBC,j,k) = 0.0_mytype
        v(iBC,j,k) = 0.0_mytype
        w(iBC,j,k) = 0.0_mytype
      enddo
    enddo
    !$acc wait(1) 
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm 
        do c=1,nHalo 
          rho( iBC+c,j,k) = 2.0_mytype*rho(iBC,j,k) - rho(iBC-c,j,k)
            u( iBC+c,j,k) = 2.0_mytype*  u(iBC,j,k) -   u(iBC-c,j,k)
            v( iBC+c,j,k) = 2.0_mytype*  v(iBC,j,k) -   v(iBC-c,j,k)
            w( iBC+c,j,k) = 2.0_mytype*  w(iBC,j,k) -   w(iBC-c,j,k)
          ien( iBC+c,j,k) = 2.0_mytype*ien(iBC,j,k) - ien(iBC-c,j,k)
          pre( iBC+c,j,k) =                           pre(iBC-c,j,k)
          tem( iBC+c,j,k) =                           tem(iBC-c,j,k)
           mu( iBC+c,j,k) = 2.0_mytype* mu(iBC,j,k) -  mu(iBC-c,j,k)
           ka( iBC+c,j,k) = 2.0_mytype* ka(iBC,j,k) -  ka(iBC-c,j,k)
        enddo
      enddo
    enddo
  ! isothermal standard boundary condition
  else if ( BC_top == "isoth_std" ) then 
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k=1,km
      do j=1,jm
        u(iBC,j,k) = 0.0_mytype
        v(iBC,j,k) = 0.0_mytype
        w(iBC,j,k) = 0.0_mytype
      enddo
    enddo
    !$acc wait(1)
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm
        pre(iBC,j,k) = 0.0_mytype
        ! prescribing the wall temperature
        tem(iBC,j,k) = Twall_top
      enddo
    enddo
    !$acc wait(1)    
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm
        !$acc loop seq
        do c = -2,-1
          pre(iBC,j,k) = pre(iBC,j,k) - c1_BD2(c)*pre(iBC+c,j,k)/c1_BD2(0)
        enddo
      enddo
    enddo
    !$acc wait(1)    
    call calcState_PT(pre,tem,rho,ien,mu,ka,iBC,iBC,1,jm,1,km)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm 
        do c=1,nHalo 
          rho( iBC+c,j,k) = 2.0_mytype*rho(iBC,j,k) - rho(iBC-c,j,k)
            u( iBC+c,j,k) = 2.0_mytype*  u(iBC,j,k) -   u(iBC-c,j,k)
            v( iBC+c,j,k) = 2.0_mytype*  v(iBC,j,k) -   v(iBC-c,j,k)
            w( iBC+c,j,k) = 2.0_mytype*  w(iBC,j,k) -   w(iBC-c,j,k)
          ien( iBC+c,j,k) = 2.0_mytype*ien(iBC,j,k) - ien(iBC-c,j,k)
          pre( iBC+c,j,k) =                           pre(iBC-c,j,k)
          tem( iBC+c,j,k) = 2.0_mytype*tem(iBC,j,k) - tem(iBC-c,j,k)
           mu( iBC+c,j,k) = 2.0_mytype* mu(iBC,j,k) -  mu(iBC-c,j,k)
           ka( iBC+c,j,k) = 2.0_mytype* ka(iBC,j,k) -  ka(iBC-c,j,k)
        enddo
      enddo
    enddo
 ! isothermal non-reflecting boundary condition
  else if (BC_top == "isoth_nrbc") then
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k=1,km
      do j=1,jm
        u(iBC,j,k) = 0.0_mytype
        v(iBC,j,k) = 0.0_mytype
        w(iBC,j,k) = 0.0_mytype
      enddo
    enddo 
    !$acc wait(1)
    !$acc parallel loop collapse(2) default(present) async(1)
    do k=1,km
      do j=1,jm 
        ! prescribing the wall temperature
        tem(iBC,j,k) = Twall_top
      enddo
    enddo
    !$acc wait(1)
    call calcState_rT(rho,tem,ien,pre,mu,ka,iBC,iBC,1,jm,1,km)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm 
        do c=1,nHalo 
          rho( iBC+c,j,k) = 2.0_mytype*rho(iBC,j,k) - rho(iBC-c,j,k)
            u( iBC+c,j,k) = 2.0_mytype*  u(iBC,j,k) -   u(iBC-c,j,k)
            v( iBC+c,j,k) = 2.0_mytype*  v(iBC,j,k) -   v(iBC-c,j,k)
            w( iBC+c,j,k) = 2.0_mytype*  w(iBC,j,k) -   w(iBC-c,j,k)
          ien( iBC+c,j,k) = 2.0_mytype*ien(iBC,j,k) - ien(iBC-c,j,k)
          pre( iBC+c,j,k) =                           pre(iBC-c,j,k)
          tem( iBC+c,j,k) = 2.0_mytype*tem(iBC,j,k) - tem(iBC-c,j,k)
           mu( iBC+c,j,k) = 2.0_mytype* mu(iBC,j,k) -  mu(iBC-c,j,k)
           ka( iBC+c,j,k) = 2.0_mytype* ka(iBC,j,k) -  ka(iBC-c,j,k)
        enddo
      enddo
    enddo
  ! non-existing boundary condition  
  else
    if (nrank.eq.0) then 
      write(stdout,*) "unknown BC_top", BC_top
    endif
    call decomp_2d_finalize
    call mpi_finalize(ierr)
    stop
  endif
end subroutine


! RHS: top boundary 
subroutine setBC_RHS_Top(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e, rho,u,v,w,ien,pre)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_eos_var
  use mod_halo
  use mod_finitediff
  implicit none
  real(mytype), dimension(:,:,:) :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre
  real(mytype), dimension(5) :: d,L
  real(mytype) :: xa, za, rhoa, ua, va, wa, iena
  real(mytype) :: Kfact, sigm, sos, fac, dp, drho, du, dv, dw, ht, prefac_r
  integer :: iBC,j,k,k1,jm,km,ierr,c
  ! at the last mesh cell in x-direction
  iBC = xsize(1)
  jm = xsize(2)
  km = xsize(3)
  k1 = 1
  ! calculate characteristics for freestream boundary condition
  if (BC_top == "free_nrbc") then 
    ! constant for the pressure prescription
    sigm = 0.25_mytype
#if defined(IG)
      prefac_r = t_ig%prefac_r 
#elif defined(VdW)
      prefac_r = t_vdw%prefac_r 
#elif defined(RK)
      prefac_r = t_rk%prefac_r
#elif defined(PR)
      prefac_r = t_pr%prefac_r 
#endif
    !$acc parallel loop collapse(2) default(present) private(sos,fac,dp,drho,du,dv,dw,d,L,Kfact,ht) async(1)
    do k = 1,km
      do j = 1,jm
        xa = xp(iBC)
        rhoa = rho(iBC,j,k)
        ua = u(iBC,j,k)
        va = v(iBC,j,k)
        wa = w(iBC,j,k)
        iena = ien(iBC,j,k)
        call calcSOS(rhoa,iena,sos)
        call calcFac(rhoa,iena,fac)
        call calc_conv_BD_ddx(dp,pre,iBC,j,k,nHalo,xa,dx)
        call calc_conv_BD_ddx(drho,rho,iBC,j,k,nHalo,xa,dx)
        call calc_conv_BD_ddx(du,u,iBC,j,k,nHalo,xa,dx)
        call calc_conv_BD_ddx(dv,v,iBC,j,k,nHalo,xa,dx)
        call calc_conv_BD_ddx(dw,w,iBC,j,k,nHalo,xa,dx)
        ! K factor according to Rudy & Strikwerda, JCP 36, 1980.
        Kfact = sigm*(1.0 - Ma**2)*sos/len_z
        ! amplitudes of characteristic waves
        L(1) = Kfact*(pre(iBC,j,k) - Pref*prefac_r)
        L(2) =  ua*(dp - (sos**2)*drho)
        L(3) =  ua*dv
        L(4) =  ua*dw
        L(5) = (ua + sos)*(dp + rhoa*sos*du)
        ! time variation of the wave amplitudes
        d(1) = (0.5*(L(5) + L(1)) - L(2))/sos**2
        d(2) = (L(5) - L(1))/(2.0*rhoa*sos)
        d(3) =  L(3)
        d(4) =  L(4)
        d(5) = fac*L(2)/sos**2
        ht = iena + pre(iBC,j,k)/rhoa + 0.5_mytype*(ua**2 + va**2 + wa**2)
        ! correct the right hand side
        rhs_r(iBC,j,k) = rhs_r(iBC,j,k) - d(1)
        rhs_u(iBC,j,k) = rhs_u(iBC,j,k) - d(1)*  ua - d(2)*rhoa
        rhs_v(iBC,j,k) = rhs_v(iBC,j,k) - d(1)*  va - d(3)*rhoa
        rhs_w(iBC,j,k) = rhs_w(iBC,j,k) - d(1)*  wa - d(4)*rhoa
        rhs_e(iBC,j,k) = rhs_e(iBC,j,k) - d(1)*ht &
                                        - d(2)*rhoa*ua &
                                        - d(3)*rhoa*va &
                                        - d(4)*rhoa*wa &
                                        - d(5) 
      enddo
    enddo
  ! calculate characteristics for adiabatic non-reflecting boundary condition
  else if (BC_top == "adiab_nrbc") then
    !$acc parallel loop collapse(2) default(present) private(sos,fac,dp,du,d,L,ht) async(1)
    do k = 1,km
      do j = 1,jm
        xa = xp(iBC)
        ua = u(iBC,j,k)
        rhoa = rho(iBC,j,k)
        iena = ien(iBC,j,k)
        call calcSOS(rhoa,iena,sos)
        call calcFac(rhoa,iena,fac)
        call calc_conv_BD_ddx(dp,pre,iBC,j,k,nHalo,xa,dx)
        call calc_conv_BD_ddx(du,u,iBC,j,k,nHalo,xa,dx)
        ! amplitudes of characteristic waves
        L(1) = (ua + sos)*(dp + rhoa*sos*du)
        L(2) = 0.0_mytype
        L(3) = 0.0_mytype
        L(4) = 0.0_mytype
        L(5) = L(1) - 2*sos*Ri_unit*(rhoa-1.0_mytype)
        ! time variation of the wave amplitudes
        d(1) = (0.5_mytype*(L(5) + L(1)) - L(2))/sos**2
        d(2) = (L(5) - L(1))/(2*rhoa*sos)
        d(3) = L(3)
        d(4) = L(4)
        d(5) = fac*L(2)/sos**2
        ht = iena + pre(iBC,j,k)/rhoa + 0.5_mytype*(ua**2 + v(iBC,j,k)**2 + w(iBC,j,k)**2)
        ! correct the right hand side
        rhs_r(iBC,j,k) = rhs_r(iBC,j,k) - d(1)
        rhs_u(iBC,j,k) = 0.0_mytype
        rhs_v(iBC,j,k) = 0.0_mytype
        rhs_w(iBC,j,k) = 0.0_mytype
        rhs_e(iBC,j,k) = rhs_e(iBC,j,k) - d(1)*ht &
                                        - d(5)
      enddo
    enddo
  ! calculate characteristics for isothermal non-reflecting boundary condition
  else if (BC_top == "isoth_nrbc") then
    !$acc parallel loop collapse(2) default(present) private(sos,fac,dp,du,d,L,ht) async(1)
    do k = 1,km
      do j = 1,jm
        xa = xp(iBC)
        ua = u(iBC,j,k)
        rhoa = rho(iBC,j,k)
        iena = ien(iBC,j,k)
        call calcSOS(rhoa,iena,sos)
        call calcFac(rhoa,iena,fac)
        call calc_conv_BD_ddx(dp,pre,iBC,j,k,nHalo,xa,dx)
        call calc_conv_BD_ddx(du,u,iBC,j,k,nHalo,xa,dx)
        ! amplitudes of characteristic waves
        L(1) = (ua + sos)*(dp + rhoa*sos*du)
        L(2) = 0.0_mytype
        L(3) = 0.0_mytype
        L(4) = 0.0_mytype
        L(5) = L(1) - 2*sos*Ri_unit*(rhoa-1.0_mytype)
        ! time variation of the wave amplitudes
        d(1) = (0.5_mytype*(L(5) + L(1)) - L(2))/sos**2
        d(2) = (L(5) - L(1))/(2*rhoa*sos)
        d(3) = L(3)
        d(4) = L(4)
        d(5) = fac*L(2)/sos**2
        ! correct the right hand side
        rhs_r(iBC,j,k) = rhs_r(iBC,j,k) - d(1)
        rhs_u(iBC,j,k) = 0.0_mytype
        rhs_v(iBC,j,k) = 0.0_mytype
        rhs_w(iBC,j,k) = 0.0_mytype
        rhs_e(iBC,j,k) = 0.0_mytype
      enddo
    enddo 
  else
    !$acc parallel loop collapse(2) default(present) async(1) 
    do k=1,km
      do j=1,jm
        rhs_r(iBC,j,k) = 0.0_mytype
        rhs_u(iBC,j,k) = 0.0_mytype
        rhs_v(iBC,j,k) = 0.0_mytype
        rhs_w(iBC,j,k) = 0.0_mytype
        rhs_e(iBC,j,k) = 0.0_mytype
      enddo
    enddo
  endif
end subroutine


! halo cells: inlet boundary
subroutine setBC_Inl(rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_halo
  use mod_finitediff 
  use mod_grid
  use mod_init
  implicit none
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl
  real(mytype), dimension(3) :: b_w = (/3.0_mytype, -3.0_mytype, 1.0_mytype/) 
  real(mytype), dimension(3) :: b_o = (/1.0_mytype,  0.0_mytype, 0.0_mytype/)
  real(mytype) :: time, rwave, iwave, sum_pre
  integer :: kBC, i,j,c, ierr, im,jm
  complex(mytype) :: zz = (0.0_mytype,1.0_mytype)
  !$acc enter data copyin(b_w,b_o,zz,alphaLST,epsilonLST)
  ! at the first mesh cell in z-direction
  im = xsize(1)
  jm = xsize(2)
  kBC = 1 
  ! non-reflecting boundary condition
  if (BC_inl == "inlet_nrbc") then 
    !$acc parallel loop collapse(3) default(present) async(1)
    do j=1,jm
      do c=1,nHalo 
        do i=1,im
          rho(i,j,kBC-c) =  2.0_mytype*rho(i,j,kBC) - rho(i,j,kBC+c)
            u(i,j,kBC-c) =  2.0_mytype*  u(i,j,kBC) -   u(i,j,kBC+c)
            v(i,j,kBC-c) =  2.0_mytype*  v(i,j,kBC) -   v(i,j,kBC+c)
            w(i,j,kBC-c) =  2.0_mytype*  w(i,j,kBC) -   w(i,j,kBC+c)
          ien(i,j,kBC-c) =  2.0_mytype*ien(i,j,kBC) - ien(i,j,kBC+c)
          pre(i,j,kBC-c) =  2.0_mytype*pre(i,j,kBC) - pre(i,j,kBC+c)
          tem(i,j,kBC-c) =  2.0_mytype*tem(i,j,kBC) - tem(i,j,kBC+c)
           mu(i,j,kBC-c) =  2.0_mytype* mu(i,j,kBC) -  mu(i,j,kBC+c)
           ka(i,j,kBC-c) =  2.0_mytype* ka(i,j,kBC) -  ka(i,j,kBC+c)
        enddo
      enddo
    enddo
  ! standard boundary condition
  else if (BC_inl == "inlet_std") then 
    !$acc parallel loop collapse(2) default(present) async(1)
    do j=1,jm
     do i=1,im 
       pre(i,j,kBC) = 0.0_mytype
     enddo
    enddo
    !$acc wait(1)
    ! pressure extrapolation: see Wasistho et al., Comput. Fluids 26, 1997
    !$acc parallel loop collapse(2) default(present) private(sum_pre) async(1)
    do j=1,jm
       do i=1,im
         sum_pre=0.0_mytype
         !$acc loop seq
         do c=1,3 
         sum_pre = sum_pre + ((1.0_mytype-w_bl(i,j,kBC))*b_w(c) + w_bl(i,j,kBC)*b_o(c))/1.0_mytype*pre(i,j,kBC+c) 
         enddo
       pre(i,j,kBC)=sum_pre
     enddo
    enddo
    !$acc wait(1)
    call calcState_PT(pre,tem,rho,ien,mu,ka,1,xsize(1),1,xsize(2),kBC,kBC)
    !$acc wait(1)
    !$acc parallel loop collapse(3) default(present) async(1)
    do j=1,jm
      do c=1,nHalo 
        do i=1,im
          rho(i,j,kBC-c) =  2.0_mytype*rho(i,j,kBC) - rho(i,j,kBC+c)
            u(i,j,kBC-c) =  2.0_mytype*  u(i,j,kBC) -   u(i,j,kBC+c)
            v(i,j,kBC-c) =  2.0_mytype*  v(i,j,kBC) -   v(i,j,kBC+c)
            w(i,j,kBC-c) =  2.0_mytype*  w(i,j,kBC) -   w(i,j,kBC+c)
          ien(i,j,kBC-c) =  2.0_mytype*ien(i,j,kBC) - ien(i,j,kBC+c)
          pre(i,j,kBC-c) =  2.0_mytype*pre(i,j,kBC) - pre(i,j,kBC+c)
          tem(i,j,kBC-c) =  2.0_mytype*tem(i,j,kBC) - tem(i,j,kBC+c)
           mu(i,j,kBC-c) =  2.0_mytype* mu(i,j,kBC) -  mu(i,j,kBC+c)
           ka(i,j,kBC-c) =  2.0_mytype* ka(i,j,kBC) -  ka(i,j,kBC+c)
        enddo
      enddo
    enddo
  ! prescribing eigenfunctions based on Linear Stability Theory
  else if (BC_inl == "inlet_lst") then 
    !$acc parallel loop collapse(2) default(present) async(1)
    do j=1,jm
      do i=1,im 
      pre(i,j,kBC) = 0.0_mytype
      w(i,j,kBC) = 0.0_mytype
      u(i,j,kBC) = 0.0_mytype
      tem(i,j,kBC) = 0.0_mytype
      enddo
    enddo
    !$acc kernels default(present) async(1)
    rwave=real(exp(zz*(alphaLST*zStartDNS-omega1*time)))
    iwave=aimag(exp(zz*(alphaLST*zStartDNS-omega1*time)))
    !$acc end kernels
    !$acc wait(1)
    ! boundary layer profiles plus Fourier ansatz, epsilonLST is an input parameter
    !$acc parallel loop collapse(2) default(present) async(1)
    do j=1,jm
      do i=1,im
        w(i,j,kBC)   = w_bl(i,j,kBC)   + epsilonLST*Ma*(w_LST_real(i,j)*rwave-w_LST_imag(i,j)*iwave)
        u(i,j,kBC)   = u_bl(i,j,kBC)   + epsilonLST*Ma*(u_LST_real(i,j)*rwave-u_LST_imag(i,j)*iwave) !  
        tem(i,j,kBC) = tem_bl(i,j,kBC) + epsilonLST*Ma*(tem_LST_real(i,j)*rwave-tem_LST_imag(i,j)*iwave)
      enddo
    enddo
    !$acc parallel loop collapse(2) default(present) private(sum_pre) async(1)
    do j=1,jm
       do i=1,im
         sum_pre=0.0_mytype
         !$acc loop seq
         do c=1,3 
         sum_pre = sum_pre + ((1.0_mytype-w_bl(i,j,kBC))*b_w(c) + w_bl(i,j,kBC)*b_o(c))/1.0_mytype*pre(i,j,kBC+c) 
         enddo
       pre(i,j,kBC)=sum_pre
     enddo
    enddo
    !$acc wait(1)
    call calcState_PT(pre,tem,rho,ien,mu,ka, 1,xsize(1),1,xsize(2),kBC,kBC)
    !$acc wait(1)
    !$acc parallel loop collapse(3) default(present) async(1)
    do j=1,jm
      do c=1,nHalo 
        do i=1,im
          rho(i,j,kBC-c) =  2.0_mytype*rho(i,j,kBC) - rho(i,j,kBC+c)
            u(i,j,kBC-c) =  2.0_mytype*  u(i,j,kBC) -   u(i,j,kBC+c)
            v(i,j,kBC-c) =  2.0_mytype*  v(i,j,kBC) -   v(i,j,kBC+c)
            w(i,j,kBC-c) =  2.0_mytype*  w(i,j,kBC) -   w(i,j,kBC+c)
          ien(i,j,kBC-c) =  2.0_mytype*ien(i,j,kBC) - ien(i,j,kBC+c)
          pre(i,j,kBC-c) =  2.0_mytype*pre(i,j,kBC) - pre(i,j,kBC+c)
          tem(i,j,kBC-c) =  2.0_mytype*tem(i,j,kBC) - tem(i,j,kBC+c)
           mu(i,j,kBC-c) =  2.0_mytype* mu(i,j,kBC) -  mu(i,j,kBC+c)
           ka(i,j,kBC-c) =  2.0_mytype* ka(i,j,kBC) -  ka(i,j,kBC+c)
        enddo
      enddo
    enddo
    ! isothermal standard boundary condition
  else if ( BC_inl == "isoth_std" ) then 
    !$acc parallel loop collapse(2) default(present) async(1) 
    do j=1,jm
      do i=1,im
        u(i,j,kBC) = 0.0_mytype
        v(i,j,kBC) = 0.0_mytype
        w(i,j,kBC) = 0.0_mytype
      enddo
    enddo
    !$acc wait(1)
    !$acc parallel loop collapse(2) default(present) async(1)
    do j=1,jm
      do i=1,im
        pre(i,j,kBC) = 0.0_mytype
        ! prescribing the inlet temperature
        tem(i,j,kBC) = Twall_inl
      enddo
    enddo
    !$acc wait(1)    
    !$acc parallel loop collapse(2) default(present) async(1)
    do j=1,jm
      do i=1,im
        !$acc loop seq
        do c = 1,2
          pre(i,j,kBC) = pre(i,j,kBC) - c1_FD2(c)*pre(i,j,kBC+c)/c1_FD2(0)
        enddo
      enddo
    enddo
    !$acc wait(1)    
    call calcState_PT(pre,tem,rho,ien,mu,ka,1,im,1,jm,kBC,kBC)
    !$acc parallel loop collapse(3) default(present) async(1)
    do j=1,jm
      do i=1,im 
        do c=1,nHalo 
          rho( i,j,kBC-c) = 2.0_mytype*rho(i,j,kBC) - rho(i,j,kBC+c)
            u( i,j,kBC-c) = 2.0_mytype*  u(i,j,kBC) -   u(i,j,kBC+c)
            v( i,j,kBC-c) = 2.0_mytype*  v(i,j,kBC) -   v(i,j,kBC+c)
            w( i,j,kBC-c) = 2.0_mytype*  w(i,j,kBC) -   w(i,j,kBC+c)
          ien( i,j,kBC-c) = 2.0_mytype*ien(i,j,kBC) - ien(i,j,kBC+c)
          pre( i,j,kBC-c) =                           pre(i,j,kBC+c)
          tem( i,j,kBC-c) = 2.0_mytype*tem(i,j,kBC) - tem(i,j,kBC+c)
           mu( i,j,kBC-c) = 2.0_mytype* mu(i,j,kBC) -  mu(i,j,kBC+c)
           ka( i,j,kBC-c) = 2.0_mytype* ka(i,j,kBC) -  ka(i,j,kBC+c)
        enddo
      enddo
    enddo
  ! non-existing boundary condition
  else
    if (nrank.eq.0) write(stdout,*) "unknown BC_inlet", BC_inl
    call decomp_2d_finalize
    call mpi_finalize(ierr)
    stop
  endif
end subroutine


! rescaling the inlet boundary condition
subroutine setBC_Inl_rescale(rho,u,v,w,ien,pre,tem,mu,ka)
  use mod_param
  use mod_eos
  use mod_halo
  use mod_grid  
  use mod_math
  implicit none
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: u1d_flu,v1d_flu,w1d_flu,p1d_flu,t1d_flu
  !$acc declare create(u1d_flu,v1d_flu,w1d_flu,p1d_flu,t1d_flu)
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: aw_rcy_vd
  !$acc declare create(aw_rcy_vd)
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: u1d_rcy_flu,v1d_rcy_flu,w1d_rcy_flu, &
                                                                                        p1d_rcy_flu,t1d_rcy_flu
  !$acc declare create(u1d_rcy_flu,v1d_rcy_flu,w1d_rcy_flu,p1d_rcy_flu,t1d_rcy_flu)
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: u1d_rcy_ave,v1d_rcy_ave,w1d_rcy_ave, &
                                                                                        p1d_rcy_ave,t1d_rcy_ave
  !$acc declare create(u1d_rcy_ave,v1d_rcy_ave,w1d_rcy_ave,p1d_rcy_ave,t1d_rcy_ave)                                                                                      
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: u_rcy_inner_flu,v_rcy_inner_flu, &
                                                                                        w_rcy_inner_flu,p_rcy_inner_flu, &
                                                                                        t_rcy_inner_flu
  !$acc declare create(u_rcy_inner_flu,v_rcy_inner_flu,w_rcy_inner_flu,p_rcy_inner_flu,t_rcy_inner_flu)                                                                                        
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: u_rcy_inner_ave,v_rcy_inner_ave, &
                                                                                        w_rcy_inner_ave,p_rcy_inner_ave, &
                                                                                        t_rcy_inner_ave
  !$acc declare create(u_rcy_inner_ave,v_rcy_inner_ave,w_rcy_inner_ave,p_rcy_inner_ave,t_rcy_inner_ave)                                                                              
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: u_rcy_outer_flu,v_rcy_outer_flu, & 
                                                                                        w_rcy_outer_flu,p_rcy_outer_flu, &
                                                                                        t_rcy_outer_flu
  !$acc declare create(u_rcy_outer_flu,v_rcy_outer_flu,w_rcy_outer_flu,p_rcy_outer_flu,t_rcy_outer_flu)                                                                                        
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: u_rcy_outer_ave,v_rcy_outer_ave, &
                                                                                        w_rcy_outer_ave,p_rcy_outer_ave, &
                                                                                        t_rcy_outer_ave
  !$acc declare create(u_rcy_outer_ave,v_rcy_outer_ave,w_rcy_outer_ave,p_rcy_outer_ave,t_rcy_outer_ave)                                                                                        
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: u_inl_inner,v_inl_inner,w_inl_inner, &
                                                                                        p_inl_inner,t_inl_inner
  !$acc declare create(u_inl_inner,v_inl_inner,w_inl_inner,p_inl_inner,t_inl_inner)  
  real(mytype), dimension(1:xsize(1), 1:xsize(2), 1-nHalo:1) :: u_inl_outer,v_inl_outer,w_inl_outer, &
                                                                                        p_inl_outer,t_inl_outer
  !$acc declare create(u_inl_outer,v_inl_outer,w_inl_outer,p_inl_outer,t_inl_outer)      
  real(mytype), dimension(1:xsize(1)) :: tmp_u_in, tmp_u_out, tmp_v_in, tmp_v_out, tmp_w_in, tmp_w_out, &
                                         tmp_p_in, tmp_p_out, tmp_t_in, tmp_t_out                                                                                 
  integer :: i,j,k
  ! sending variables from rescaling location to inlet
  if ((neigh%inlet) .or. (neigh%recycle)) then
     call haloUpdate_rescale(xsize,u,v,w,pre,tem)
  endif
  ! rescaling the inlet boundary layer thickness
  if (neigh%inlet) then
    !$acc parallel loop collapse(3) default(present)
    do k=1-nHalo,1
      do j=1,xsize(2)
        do i=1,xsize(1)
          ! fluctuations
          u1d_flu(i,j,k) = u(i,j,k)   - au_rcy(i)
          v1d_flu(i,j,k) = v(i,j,k)   - av_rcy(i)
          w1d_flu(i,j,k) = w(i,j,k)   - aw_rcy(i)
          p1d_flu(i,j,k) = pre(i,j,k) - ap_rcy(i)
          t1d_flu(i,j,k) = tem(i,j,k) - at_rcy(i)
          ! averaging (Van-Driest transformation for streamwise velocity with u_inf as 1)
          aw_rcy_vd(i,j,k) = 1.0_mytype/A_FF * asin(A_FF*aw_rcy(i)) 
        enddo
      enddo
    enddo
    ! interpolation fluctuation
    !$acc parallel loop collapse(2) default(present) &
    !$acc private(tmp_u_in,tmp_v_in,tmp_w_in,tmp_p_in,tmp_t_in) &
    !$acc private(tmp_u_out,tmp_v_out,tmp_w_out,tmp_p_out,tmp_t_out)
    do k=1-nHalo,1
      do j=1,xsize(2)
        tmp_u_in = u1d_flu(1:xsize(1),j,k)
        tmp_v_in = v1d_flu(1:xsize(1),j,k)
        tmp_w_in = w1d_flu(1:xsize(1),j,k)
        tmp_p_in = p1d_flu(1:xsize(1),j,k)
        tmp_t_in = t1d_flu(1:xsize(1),j,k)
        call spline(x,tmp_u_in,xsize(1),tmp_u_out)
        call spline(x,tmp_v_in,xsize(1),tmp_v_out)
        call spline(x,tmp_w_in,xsize(1),tmp_w_out)
        call spline(x,tmp_p_in,xsize(1),tmp_p_out)
        call spline(x,tmp_t_in,xsize(1),tmp_t_out)
        u1d_rcy_flu(1:xsize(1),j,k) = tmp_u_out
        v1d_rcy_flu(1:xsize(1),j,k) = tmp_v_out
        w1d_rcy_flu(1:xsize(1),j,k) = tmp_w_out
        p1d_rcy_flu(1:xsize(1),j,k) = tmp_p_out
        t1d_rcy_flu(1:xsize(1),j,k) = tmp_t_out
      enddo
    enddo
    ! interpolation average
    !$acc parallel loop collapse(2) default(present) &
    !$acc private(tmp_w_in) &
    !$acc private(tmp_u_out,tmp_v_out,tmp_w_out,tmp_p_out,tmp_t_out)
    do k=1-nHalo,1
      do j=1,xsize(2)
        tmp_w_in = aw_rcy_vd(1:xsize(1),j,k)
        call spline(x,au_rcy   ,xsize(1),tmp_u_out)
        call spline(x,av_rcy   ,xsize(1),tmp_v_out)
        call spline(x,tmp_w_in ,xsize(1),tmp_w_out)
        call spline(x,ap_rcy   ,xsize(1),tmp_p_out)
        call spline(x,at_rcy   ,xsize(1),tmp_t_out)
        u1d_rcy_ave(1:xsize(1),j,k) = tmp_u_out
        v1d_rcy_ave(1:xsize(1),j,k) = tmp_v_out
        w1d_rcy_ave(1:xsize(1),j,k) = tmp_w_out
        p1d_rcy_ave(1:xsize(1),j,k) = tmp_p_out
        t1d_rcy_ave(1:xsize(1),j,k) = tmp_t_out
      enddo
    enddo
    !$acc parallel loop collapse(3) default(present)
    do k=1-nHalo,1
      do j=1,xsize(2)
        do i=1,xsize(1)
          ! interpolating the inner layer
          ! fluctuations
          call splint(x,u1d_flu(1:xsize(1),j,k),u1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),u_rcy_inner_flu(i,j,k))
          call splint(x,v1d_flu(1:xsize(1),j,k),v1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),v_rcy_inner_flu(i,j,k))
          call splint(x,w1d_flu(1:xsize(1),j,k),w1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),w_rcy_inner_flu(i,j,k))
          call splint(x,p1d_flu(1:xsize(1),j,k),p1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),p_rcy_inner_flu(i,j,k))
          call splint(x,t1d_flu(1:xsize(1),j,k),t1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),t_rcy_inner_flu(i,j,k))
          ! fluctuations are zero in the free stream
          if (x(i) .gt. delta_inl*1.5_mytype) then
            u_rcy_inner_flu(i,j,k) = 0.0_mytype
            v_rcy_inner_flu(i,j,k) = 0.0_mytype
            w_rcy_inner_flu(i,j,k) = 0.0_mytype
            p_rcy_inner_flu(i,j,k) = 0.0_mytype
            t_rcy_inner_flu(i,j,k) = 0.0_mytype
          endif
          ! average
          call splint(x,au_rcy(1:xsize(1))   ,u1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),u_rcy_inner_ave(i,j,k))
          call splint(x,av_rcy(1:xsize(1))   ,v1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),v_rcy_inner_ave(i,j,k))
          call splint(x,aw_rcy_vd(1:xsize(1),j,k),w1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),w_rcy_inner_ave(i,j,k))
          call splint(x,ap_rcy(1:xsize(1))   ,p1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),p_rcy_inner_ave(i,j,k))
          call splint(x,at_rcy(1:xsize(1))   ,t1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_inner(i),t_rcy_inner_ave(i,j,k))
          ! use free stream values
          if (x(i) .gt. delta_inl*1.5_mytype) then
            u_rcy_inner_ave(i,j,k) = 0.0_mytype 
            v_rcy_inner_ave(i,j,k) = 0.0_mytype 
            w_rcy_inner_ave(i,j,k) = 1.0_mytype 
            p_rcy_inner_ave(i,j,k) = Pref       
            t_rcy_inner_ave(i,j,k) = 1.0_mytype
          endif    
          ! interpolating the outer layer
          ! fluctuations
          call splint(x,u1d_flu(1:xsize(1),j,k),u1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),u_rcy_outer_flu(i,j,k))
          call splint(x,v1d_flu(1:xsize(1),j,k),v1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),v_rcy_outer_flu(i,j,k))
          call splint(x,w1d_flu(1:xsize(1),j,k),w1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),w_rcy_outer_flu(i,j,k))
          call splint(x,p1d_flu(1:xsize(1),j,k),p1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),p_rcy_outer_flu(i,j,k))
          call splint(x,t1d_flu(1:xsize(1),j,k),t1d_rcy_flu(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),t_rcy_outer_flu(i,j,k))
          ! fluctuations are zero in the free stream
          if (x(i) .gt. delta_inl*1.5_mytype) then
            u_rcy_outer_flu(i,j,k) = 0.0_mytype
            v_rcy_outer_flu(i,j,k) = 0.0_mytype
            w_rcy_outer_flu(i,j,k) = 0.0_mytype
            p_rcy_outer_flu(i,j,k) = 0.0_mytype
            t_rcy_outer_flu(i,j,k) = 0.0_mytype
          endif
          ! average
          call splint(x,au_rcy(1:xsize(1))   ,u1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),u_rcy_outer_ave(i,j,k))
          call splint(x,av_rcy(1:xsize(1))   ,v1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),v_rcy_outer_ave(i,j,k))
          call splint(x,aw_rcy_vd(1:xsize(1),j,k),w1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),w_rcy_outer_ave(i,j,k))
          call splint(x,ap_rcy(1:xsize(1))   ,p1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),p_rcy_outer_ave(i,j,k))
          call splint(x,at_rcy(1:xsize(1))   ,t1d_rcy_ave(1:xsize(1),j,k),xsize(1),x_rcy_outer(i),t_rcy_outer_ave(i,j,k))
          ! use free stream values
          if (x(i) .gt. delta_inl*1.5_mytype) then
            u_rcy_outer_ave(i,j,k) = 0.0_mytype 
            v_rcy_outer_ave(i,j,k) = 0.0_mytype 
            w_rcy_outer_ave(i,j,k) = 1.0_mytype 
            p_rcy_outer_ave(i,j,k) = Pref       
            t_rcy_outer_ave(i,j,k) = 1.0_mytype 
          endif
          ! addition inner averages and fluctuations
          u_inl_inner(i,j,k) = u_rcy_inner_ave(i,j,k) + beta_rescale*u_rcy_inner_flu(i,j,k)
          v_inl_inner(i,j,k) = v_rcy_inner_ave(i,j,k) + beta_rescale*v_rcy_inner_flu(i,j,k)
          w_inl_inner(i,j,k) = beta_rescale*(w_rcy_inner_ave(i,j,k) + w_rcy_inner_flu(i,j,k))
          p_inl_inner(i,j,k) = p_rcy_inner_ave(i,j,k) + p_rcy_inner_flu(i,j,k)
          t_inl_inner(i,j,k) = t_rcy_inner_ave(i,j,k) + t_rcy_inner_flu(i,j,k)
          ! addition outer averages and fluctuations
          u_inl_outer(i,j,k) = u_rcy_outer_ave(i,j,k) + beta_rescale*u_rcy_outer_flu(i,j,k)
          v_inl_outer(i,j,k) = v_rcy_outer_ave(i,j,k) + beta_rescale*v_rcy_outer_flu(i,j,k)
          w_inl_outer(i,j,k) = aw_rcy_vd(xsize(1),j,k) -beta_rescale* &
                               (aw_rcy_vd(xsize(1),j,k) - w_rcy_outer_ave(i,j,k)) + beta_rescale*w_rcy_outer_flu(i,j,k)
          p_inl_outer(i,j,k) = p_rcy_outer_ave(i,j,k) + p_rcy_outer_flu(i,j,k)
          t_inl_outer(i,j,k) = t_rcy_outer_ave(i,j,k) + t_rcy_outer_flu(i,j,k)
          ! weighting for total profiles
          u(i,j,k) = u_inl_inner(i,j,k) * (1.0_mytype - weight_func(i)) + &
                     u_inl_outer(i,j,k) * weight_func(i)
          v(i,j,k) = v_inl_inner(i,j,k) * (1.0_mytype - weight_func(i)) + &
                     v_inl_outer(i,j,k) * weight_func(i)
          w(i,j,k) = w_inl_inner(i,j,k) * (1.0_mytype - weight_func(i)) + &
                     w_inl_outer(i,j,k) * weight_func(i)
          pre(i,j,k) = p_inl_inner(i,j,k) * (1.0_mytype - weight_func(i)) + &
                       p_inl_outer(i,j,k) * weight_func(i)
          tem(i,j,k) = t_inl_inner(i,j,k) * (1.0_mytype - weight_func(i)) + &
                       t_inl_outer(i,j,k) * weight_func(i)  
          enddo
        enddo
      enddo
     ! updating secondary non-conservative variables
     call calcState_PT(pre,tem,rho,ien,mu,ka,1,xsize(1),1,xsize(2),1-nHalo,1)
  endif
end subroutine


! RHS: inlet boundary
subroutine setBC_RHS_Inl(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e, rho,u,v,w,ien,pre)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_halo
  use mod_finitediff 
  implicit none
  integer :: kBC, i,j,c, ierr, im,jm
  real(mytype), dimension(:,:,:) :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre
  real(mytype), dimension(5) :: d,L
  real(mytype) :: sos, fac, dp, dw, ht
  real(mytype) :: za, rhoa, wa, iena
  ! at the first mesh cell in z-direction
  kBC = 1
  im  = xsize(1)
  jm  = xsize(2)
  ! non-reflecting boundary condition
  if (BC_inl == "inlet_nrbc") then 
    !$acc parallel loop collapse(2) default(present) private(sos,dp,dw,d,L,ht) async(1)
    do j = 1,jm
      ! i=1 is taken by the bottom boundary condition
      do i = 2,im  
        za   = zp(kBC)
        rhoa = rho(i,j,kBC)
        wa   = w(i,j,kBC)
        call calcSOS(rhoa,ien(i,j,kBC),sos)
        call calc_conv_FD_ddz(dp,pre,i,j,kBC,nHalo,za,dz)
        call calc_conv_FD_ddz(dw,w,i,j,kBC,nHalo,za,dz)
        ! amplitude of characteristic waves
        L(1) = (wa - sos)*(dp - rhoa*sos*dw)
        ! time variation of the wave amplitude
        d(1) = L(1)/sos**2
        ht = ien(i,j,kBC) + pre(i,j,kBC)/rhoa + 0.5_mytype*(u(i,j,kBC)**2 + v(i,j,kBC)**2 + w(i,j,kBC)**2)
        ! correct the right hand side
        rhs_r(i,j,kBC) = 0.0_mytype
        rhs_u(i,j,kBC) = 0.0_mytype
        rhs_v(i,j,kBC) = 0.0_mytype
        rhs_w(i,j,kBC) = 0.0_mytype 
        rhs_e(i,j,kBC) = rhs_e(i,j,kBC) - d(1)*ht
      enddo
    enddo
  ! standard boundary condition
  else
    !$acc parallel loop collapse(2) default(present) async(1) 
    do j=1,jm
      do i=2,im
        rhs_r(i,j,kBC) = 0.0_mytype
        rhs_u(i,j,kBC) = 0.0_mytype
        rhs_v(i,j,kBC) = 0.0_mytype
        rhs_w(i,j,kBC) = 0.0_mytype
        rhs_e(i,j,kBC) = 0.0_mytype
      enddo
    enddo
  ! if BC_inl == 'inlet_lst', no RHS is needed 
  endif
end subroutine


! halo cells: outlet boundary
subroutine setBC_Out(rho,u,v,w,ien,pre,tem,mu,ka)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_halo
  use mod_finitediff 
  implicit none
  integer :: kBC,i,j,c, ierr, im,jm,km
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
  ! at the last mesh cell in z-direction
  im  = xsize(1)
  jm  = xsize(2)
  kBC = xsize(3)
  ! non-reflecting boundary condition
  if (BC_out == "outlet_nrbc") then 
    !$acc parallel loop collapse(3) default(present) async(1)
    do j=1,jm
      do c=1,nHalo 
        do i=1,im
          rho(i,j,kBC+c) = 2.0_mytype*rho(i,j,kBC) - rho(i,j,kBC-c)
            u(i,j,kBC+c) = 2.0_mytype*  u(i,j,kBC) -   u(i,j,kBC-c)
            v(i,j,kBC+c) = 2.0_mytype*  v(i,j,kBC) -   v(i,j,kBC-c)
            w(i,j,kBC+c) = 2.0_mytype*  w(i,j,kBC) -   w(i,j,kBC-c)
          ien(i,j,kBC+c) = 2.0_mytype*ien(i,j,kBC) - ien(i,j,kBC-c)
          pre(i,j,kBC+c) = 2.0_mytype*pre(i,j,kBC) - pre(i,j,kBC-c)
          tem(i,j,kBC+c) = 2.0_mytype*tem(i,j,kBC) - tem(i,j,kBC-c)
           mu(i,j,kBC+c) = 2.0_mytype* mu(i,j,kBC) -  mu(i,j,kBC-c)
           ka(i,j,kBC+c) = 2.0_mytype* ka(i,j,kBC) -  ka(i,j,kBC-c)
        enddo
      enddo
    enddo
    ! isothermal standard boundary condition
  else if ( BC_out == "isoth_std" ) then 
    !$acc parallel loop collapse(2) default(present) async(1) 
    do j=1,jm
      do i=1,im
        u(i,j,kBC) = 0.0_mytype
        v(i,j,kBC) = 0.0_mytype
        w(i,j,kBC) = 0.0_mytype
      enddo
    enddo
    !$acc wait(1)
    !$acc parallel loop collapse(2) default(present) async(1)
    do j=1,jm
      do i=1,im
        pre(i,j,kBC) = 0.0_mytype
        ! prescribing the outlet temperature
        tem(i,j,kBC) = Twall_out
      enddo
    enddo
    !$acc wait(1)    
    !$acc parallel loop collapse(2) default(present) async(1)
    do j=1,jm
      do i=1,im
        !$acc loop seq
        do c = -2,-1
          pre(i,j,kBC) = pre(i,j,kBC) - c1_BD2(c)*pre(i,j,kBC+c)/c1_BD2(0)
        enddo
      enddo
    enddo
    !$acc wait(1)    
    call calcState_PT(pre,tem,rho,ien,mu,ka,1,im,1,jm,kBC,kBC)
    !$acc parallel loop collapse(3) default(present) async(1)
    do j=1,jm
      do i=1,im 
        do c=1,nHalo 
          rho( i,j,kBC+c) = 2.0_mytype*rho(i,j,kBC) - rho(i,j,kBC-c)
            u( i,j,kBC+c) = 2.0_mytype*  u(i,j,kBC) -   u(i,j,kBC-c)
            v( i,j,kBC+c) = 2.0_mytype*  v(i,j,kBC) -   v(i,j,kBC-c)
            w( i,j,kBC+c) = 2.0_mytype*  w(i,j,kBC) -   w(i,j,kBC-c)
          ien( i,j,kBC+c) = 2.0_mytype*ien(i,j,kBC) - ien(i,j,kBC-c)
          pre( i,j,kBC+c) =                           pre(i,j,kBC-c)
          tem( i,j,kBC+c) = 2.0_mytype*tem(i,j,kBC) - tem(i,j,kBC-c)
           mu( i,j,kBC+c) = 2.0_mytype* mu(i,j,kBC) -  mu(i,j,kBC-c)
           ka( i,j,kBC+c) = 2.0_mytype* ka(i,j,kBC) -  ka(i,j,kBC-c)
        enddo
      enddo
    enddo
  ! non-existing boundary condition
  else
      if (nrank.eq.0) write(stdout,*) "unknown BC_outlet"
      call decomp_2d_finalize
      call mpi_finalize(ierr)
  endif
end subroutine


! RHS: outlet boundary
subroutine setBC_RHS_Out(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,rho,u,v,w,ien,pre,p_ref)
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_eos_var
  use mod_halo
  use mod_finitediff 
  implicit none
  integer :: kBC, i,j,c,ierr,im,jm
  real(mytype), dimension(:,:,:) :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e  
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre
  real(mytype), dimension(5) :: d,L
  real(mytype), dimension(:,:) :: p_ref
  real(mytype) :: sos, fac, Kfact, dp, drho, du, dv, dw, ht, prefac_r
  real(mytype) :: za, rhoa, ua, va, wa, iena
  im  = xsize(1)
  jm  = xsize(2)
  kBC = xsize(3)
  ! non-reflecting boundary condition
  if (BC_out == "outlet_nrbc") then 
#if defined(IG)
    prefac_r = t_ig%prefac_r 
#elif defined(VdW)
    prefac_r = t_vdw%prefac_r
#elif defined(RK)
    prefac_r = t_rk%prefac_r 
#elif defined(PR)
    prefac_r = t_pr%prefac_r 
#endif
    !$acc parallel loop collapse(2) default(present) private(sos,fac,dp,drho,du,dv,dw,ht,Kfact,d,L,ht) async(1)
    do j = 1,jm
      ! i=1 is taken by the bottom boundary condition
      do i = 2,im 
        za = zp(kBC)
        rhoa = rho(i,j,kBC)
        ua = u(i,j,kBC)
        va = v(i,j,kBC)
        wa = w(i,j,kBC)
        iena = ien(i,j,kBC)
        call calcSOS(rhoa,iena,sos)
        call calcFac(rhoa,iena,fac)
        call calc_conv_BD_ddz(dp,pre,i,j,kBC,nHalo,za,dz)
        call calc_conv_BD_ddz(drho,rho,i,j,kBC,nHalo,za,dz)
        call calc_conv_BD_ddz(du,u,i,j,kBC,nHalo,za,dz)
        call calc_conv_BD_ddz(dv,v,i,j,kBC,nHalo,za,dz)
        call calc_conv_BD_ddz(dw,w,i,j,kBC,nHalo,za,dz)
        ! K factor according to Rudy & Strikwerda, JCP 36, 1980.
        Kfact = sigm_outlet*(1.0_mytype - Ma**2)*sos/len_z
        ! amplitudes of characteristic waves
        L(1) = Kfact*(pre(i,j,kBC) - p_ref(i,kBC)) + flag_supersonic_outlet*(wa - sos)*(dp - rhoa*sos*dw)
        L(2) =  wa*(dp - (sos**2)*drho)
        L(3) =  wa*du
        L(4) =  wa*dv
        L(5) = (wa + sos)*(dp + rhoa*sos*dw)
        ! time variation of the wave amplitude
        d(1) = (0.5_mytype*(L(5) + L(1)) - L(2))/sos**2
        d(2) = (L(5) - L(1))/(2*rhoa*sos)
        d(3) = L(3)
        d(4) = L(4)
        d(5) = fac*L(2)/sos**2
        ht = iena + pre(i,j,kBC)/rhoa + 0.5_mytype*(ua**2 + va**2 + wa**2)
        ! correcting the right hand side
        rhs_r(i,j,kBC) = rhs_r(i,j,kBC) - d(1)
        rhs_u(i,j,kBC) = rhs_u(i,j,kBC) - d(1)*ua - d(3)*rhoa
        rhs_v(i,j,kBC) = rhs_v(i,j,kBC) - d(1)*va - d(4)*rhoa
        rhs_w(i,j,kBC) = rhs_w(i,j,kBC) - d(1)*wa - d(2)*rhoa
        rhs_e(i,j,kBC) = rhs_e(i,j,kBC) - d(1)*ht &
                                        - d(2)*rhoa*wa &
                                        - d(3)*rhoa*ua &
                                        - d(4)*rhoa*va &
                                        - d(5)
      enddo
    enddo
  else
    !$acc parallel loop collapse(2) default(present) async(1) 
    do j=1,jm
      do i=2,im
        rhs_r(i,j,kBC) = 0.0_mytype
        rhs_u(i,j,kBC) = 0.0_mytype
        rhs_v(i,j,kBC) = 0.0_mytype
        rhs_w(i,j,kBC) = 0.0_mytype
        rhs_e(i,j,kBC) = 0.0_mytype
      enddo
    enddo
  endif
end subroutine


end module mod_boundary
