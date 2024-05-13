! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
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
    ! print the boundary conditions
    if (nrank == 0) then
      write(stdout,* ) 
      write(stdout,* ) 'Boundary conditions'
      write(stdout,* ) '-------------------------'
    select case (CASE)
      case("BoundaryLayer","Channel")
        if (BC_bot(1:5)==wall_bc) then
          write(stdout,* ) 'Correct BC_bot in initBL and config.h'    
        else
          write(stdout,* ) 'Mismatch BC_bot between initBL and config.h, check both again!'
          call decomp_2d_finalize
          call mpi_finalize(ierr) 
          stop
        endif
    end select
    write(stdout,* ) 'Top boundary                                         ',BC_top
    write(stdout,* ) 'Bottom boundary                                      ',BC_bot
    if (CASE == "BoundaryLayer") then
      write(stdout,* ) 'T_wall (T/T_inf)                                     ',Twall_bottom
    else if (CASE == "Channel") then
      write(stdout,* ) 'top T_wall (T/T_inf)                                     ',Twall_top
      write(stdout,* ) 'bottom T_wall (T/T_inf)                                  ',Twall_bottom
    endif
    write(stdout,* ) 'Inlet boundary                                       ',BC_inl
    if (BC_inl == "inlet_lst") then
      write(stdout,* ) '      perturbation wave at the inlet!                                '
      write(stdout,* ) 
      write(stdout,* ) '      streamwise wavenumber                             ',alphaLST
      write(stdout,* ) '      frequency                                         ',omega1
      write(stdout,* ) '      normalised amplitude                              ',epsilonLST*Ma
      write(stdout,* )
    endif
      write(stdout,* ) 'Outlet boundary                                      ',BC_out
      write(stdout,* ) '-------------------------'
      write(stdout,* ) 'Boundary calculation                                 completed!'
      write(stdout,* ) 
    endif
    if ((p_row == 1) .and. (nrank == i)) then
      if ((neigh%inlet)) then 
         write(stdout,* ) 'Inlet boundary, proc=                                       ',nrank
      endif
    endif
    if ((p_row == 1) .and. (nrank==p_col-1)) then
      if ((neigh%outlet))  then 
         write(stdout,* ) 'Outlet boundary, proc=                                      ',nrank
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
! set all boundary conditions
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
    ! top boundary
    call setBC_Bot(rho,u,v,w,ien,pre,tem,mu,ka,time)
    ! bottom boundary
    call setBC_Top(rho,u,v,w,ien,pre,tem,mu,ka)
  endif
  if ((perBC(3) .eqv. .false.) .and. (neigh%inlet)) then 
    ! inlet boundary
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
  integer :: kk
  real(mytype) :: factz, ztem_1, ztem_2
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
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm
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
        tem(iBC,j,k) = Twall_bottom
      enddo
    enddo
    !$acc wait(1)    
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm
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
        tem(iBC,j,k) = Twall_bottom
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
    if (nrank.eq.0) write (*,*) "unknown BC_bottom: ", BC_bot
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
  integer :: iBC,j,k,jm,km,c, ierr, k1, kk
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
              L(5) = L(1) - 2*rhoa*sos*dudt_pert(j,k)
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
        L(5) = L(1)
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
          L(5) = L(1) - 2*rhoa*sos*dudt_pert(j,k)
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
          L(5) = L(1)
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
  ! without boundary condition
  else 
    if (nrank.eq.0) write (*,*) "unknown BC_bottom"
    call decomp_2d_finalize
    call mpi_finalize(ierr)
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
        dudt_pert(j,k) = 0.0_mytype 
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
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm
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
        dudt_pert(j,k) = 0.0_mytype 
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
        dudt_pert(j,k) = 0.0_mytype 
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
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm
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
        dudt_pert(j,k) = 0.0_mytype 
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
    if (nrank.eq.0) write (*,*) "unknown BC_top"
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
  real(mytype) :: xa, rhoa, ua, va, wa, iena
  real(mytype) :: Kfact, sigm, sos, fac, dp, drho, du, dv, dw, ht, prefac_r
  integer :: iBC,j,k,jm,km,ierr,c
  ! at the last mesh cell in x-direction
  iBC = xsize(1)
  jm = xsize(2)
  km = xsize(3)
  ! calculate characteristics for freestream boundary condition
  if (BC_top == "free_nrbc") then 
    ! constant for the pressure prescription
    sigm = 0.25_mytype
    select case (t_param%USE_EOS)
    case("IG")
      prefac_r = t_ig%prefac_r 
    case("VdW")
      prefac_r = t_vdw%prefac_r 
    case("PR")
      prefac_r = t_pr%prefac_r 
    end select
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
        ht = iena + pre(iBC,j,k)/rhoa + 0.5*(ua**2 + va**2 + wa**2)
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
        L(1) = (ua - sos)*(dp - rhoa*sos*du)
        L(2) = 0.0_mytype
        L(3) = 0.0_mytype
        L(4) = 0.0_mytype
        L(5) = L(1)
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
        rhoa = rho(iBC,j,k)
        iena = ien(iBC,j,k)
        call calcSOS(rhoa,iena,sos)
        call calcFac(rhoa,iena,fac)
        call calc_conv_BD_ddx(dp,pre,iBC,j,k,nHalo,xa,dx)
        call calc_conv_BD_ddx(du,u,iBC,j,k,nHalo,xa,dx)
        ! amplitudes of characteristic waves
        L(1) = (u(iBC,j,k) - sos)*(dp - rhoa*sos*du)
        L(2) = 0.0_mytype
        L(3) = 0.0_mytype
        L(4) = 0.0_mytype
        L(5) = L(1)
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
  ! without boundary condition
  else 
    if (nrank.eq.0) write (*,*) "unknown BC_top"
    call decomp_2d_finalize
    call mpi_finalize(ierr)
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
    !$acc parallel loop collapse(3) default(present) 
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
  ! non-existing boundary condition
  else
    if (nrank.eq.0) write (*,*) "unknown BC_inlet"
    call decomp_2d_finalize
    call mpi_finalize(ierr)
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
  real(mytype) :: d,L
  real(mytype) :: sos, fac, dp, dw, ht
  real(mytype) :: za, rhoa
  ! at the first mesh cell in z-direction
  kBC = 1
  im  = xsize(1)
  jm  = xsize(2)
  ! non-reflecting boundary condition
  if (BC_inl == "inlet_nrbc") then 
    !$acc parallel loop collapse(2) default(present) private(sos,dp,dw,d,L,ht) async(1)
    do j = 1,xsize(2)
      ! i=1 is taken by the bottom boundary condition
      do i = 2,xsize(1)  
        za = zp(kBC)
        rhoa = rho(i,j,kBC)
        call calcSOS(rhoa,ien(i,j,kBC),sos)
        call calc_conv_FD_ddz(dp,pre,i,j,kBC,nHalo,za,dz)
        call calc_conv_FD_ddz(dw,w,i,j,kBC,nHalo,za,dz)
        ! amplitude of characteristic waves
        L = (w(i,j,kBC) - sos)*(dp - rhoa*sos*dw)
        ! time variation of the wave amplitude
        d = L/sos**2
        ht = ien(i,j,kBC) + pre(i,j,kBC)/rhoa + 0.5_mytype*(u(i,j,kBC)**2 + v(i,j,kBC)**2 + w(i,j,kBC)**2)
        ! correct the right hand side
        rhs_r(i,j,kBC) = 0.0_mytype
        rhs_u(i,j,kBC) = 0.0_mytype
        rhs_v(i,j,kBC) = 0.0_mytype
        rhs_w(i,j,kBC) = 0.0_mytype 
        rhs_e(i,j,kBC) = rhs_e(i,j,kBC) - d*ht
      enddo
    enddo
  ! standard boundary condition
  else if (BC_inl == "inlet_std") then  
    !$acc parallel loop collapse(2) default(present) async(1)
    do j=1,jm
      do i=1,im
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
  ! non-existing boundary condition
  else
      if (nrank.eq.0) write (*,*) "unknown BC_outlet"
      call decomp_2d_finalize
      call mpi_finalize(ierr)
  endif

end subroutine
! RHS: outlet boundary
subroutine setBC_RHS_Out(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,rho,u,v,w,ien,pre)
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
  real(mytype) :: sos, fac, Kfact, dp, drho, du, dv, dw, ht, prefac_r
  real(mytype) :: za, rhoa, ua, va, wa, iena
  im  = xsize(1)
  jm  = xsize(2)
  kBC = xsize(3)
  ! non-reflecting boundary condition
  if (BC_inl == "inlet_nrbc") then 
    select case (t_param%USE_EOS)
      case("IG")
        prefac_r = t_ig%prefac_r 
      case("VdW")
        prefac_r = t_vdw%prefac_r 
      case("PR")
        prefac_r = t_pr%prefac_r 
    end select
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
        Kfact = sigm_outlet*(1.0 - Ma**2)*sos/len_z
        ! amplitudes of characteristic waves
        L(1) = Kfact*(pre(i,j,kBC) - Pref*prefac_r) + flag_supersonic_outlet*(wa - sos)*(dp - rhoa*sos*dw)
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
  ! without boundary condition
  else 
    if (nrank.eq.0) write (*,*) "unknown BC_bottom"
    call decomp_2d_finalize
    call mpi_finalize(ierr)
  endif
end subroutine
end module mod_boundary
