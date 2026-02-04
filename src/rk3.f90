! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! Explicit TVD 3rd-order Runge Kutta in the low-storage version

module mod_solve
  use decomp_2d
  implicit none
  real(mytype), allocatable, dimension(:,:,:) :: drho,drhu,drhv,drhw,dret, &
                                                 rhoOld,rhuOld,rhvOld,rhwOld,retOld
contains

! initialization low-storage 3rd-order Runge-Kutta
  subroutine init_rk3()
    use decomp_2d
    use mod_param
    implicit none
    ! RHS
    allocate( drho(xsize(1), xsize(2), xsize(3)) )
    allocate( drhu(xsize(1), xsize(2), xsize(3)) )
    allocate( drhv(xsize(1), xsize(2), xsize(3)) )
    allocate( drhw(xsize(1), xsize(2), xsize(3)) )
    allocate( dret(xsize(1), xsize(2), xsize(3)) )
    ! old timestep
    allocate(rhoOld(xsize(1), xsize(2), xsize(3)) )
    allocate(rhuOld(xsize(1), xsize(2), xsize(3)) )
    allocate(rhvOld(xsize(1), xsize(2), xsize(3)) )
    allocate(rhwOld(xsize(1), xsize(2), xsize(3)) )
    allocate(retOld(xsize(1), xsize(2), xsize(3)) )
    !$acc enter data create(drho,drhu,drhv,drhw,dret)
    !$acc enter data create(rhoOld,rhuOld,rhvOld,rhwOld,retOld)
  end subroutine


! calculation of the timestep
  subroutine calcTimeStep(dt,CFL_new,rho,u,v,w,ien,mu,ka)
    use decomp_2d
    use mpi
    use mod_eos
    use mod_param
    use mod_grid
    use mod_perturbation
    implicit none
    integer :: i,j,k,ierr
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,mu,ka
    real(mytype) :: dt, dtConvInv, dtViscInv, dtCondInv, tmp, sos, delx, delz, CFL_new
    dtConvInv = 0.0_mytype
    dtViscInv = 0.0_mytype
    dtCondInv = 0.0_mytype
    ! calculation of the three (convection, viscous, conduction) timesteps
    !$acc parallel loop collapse(3) default(present) private(delx,delz,sos) reduction(max:dtConvInv,dtViscInv,dtCondInv)
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          delx = dx/xp(i)
          delz = dz/zp(k)
          call calcSOS(rho(i,j,k),ien(i,j,k),sos)
          dtConvInv =  max(dtConvInv, &
                      (abs(u(i,j,k)) + sos)/delx, &
                      (abs(v(i,j,k)) + sos)/dy,   &
                      (abs(w(i,j,k)) + sos)/delz )
          dtViscInv =  max(dtViscInv, mu(i,j,k)/delx**2,  &
                                      mu(i,j,k)/dy**2,    &
                                      mu(i,j,k)/delz**2)
          dtCondInv =  max(dtCondInv, ka(i,j,k)/delx**2,  &
                                      ka(i,j,k)/dy**2,    &
                                      ka(i,j,k)/delz**2)
        enddo
      enddo
    enddo
    !$acc wait
    ! timestep without perturbation, CFL is an input parameter
    if (pert_calc==0) then
      dt = CFL/max(dtConvInv,dtViscInv,dtCondInv)
      call mpi_allreduce(MPI_IN_PLACE, dt, 1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
      if (dtMax > 0) then
        dt = min(dt,dtMax)
      endif
    endif
    ! timestep with perturbation, a new CFL is calculated from dt_FFT (perturbation)
    if ((pert_calc==1) .OR. (BC_inl == "inlet_lst")) then
      CFL_new=dt_FFT * max(dtConvInv, dtViscInv, dtCondInv)
      call mpi_allreduce(MPI_IN_PLACE, CFL_new, 1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    else
      CFL_new=CFL
    endif
    ! the new CFL should not be larger than the input CFL 
    if (CFL_new .gt. CFL) then
      if  (nrank .eq. 0) then
        write(*,*) "CFL_FFT:", CFL_new
        write(*,*) "ERROR! Time step too big!"
      endif
      call decomp_2d_finalize
      call mpi_finalize(ierr) 
      stop
    endif
    ! a constant timestep is set when the flow is pertubed
    if ((pert_calc==1) .OR. (BC_inl == "inlet_lst")) then
      dt=dt_FFT
    endif
  end subroutine


! calculation of the three Runge-Kutta substeps and next timestep
  subroutine rk3(part,dt,istep,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
    use decomp_2d
    use mod_param
    use mod_grid
    use mod_eos
    use mod_halo
    use mod_rhs
    use mod_boundary
    implicit none
    integer istep,i,j,k,im,jm,km, ierr
    real(mytype) :: dt,time,onefrth,onethrd,twothrd,thrfrth
    real(mytype) :: rho1,u1,v1,w1,ien1,rho2,u2,v2,w2,ien2
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl
    TYPE (DECOMP_INFO), intent(IN) :: part

    onefrth = 0.25_mytype
    onethrd = 1.0_mytype/3.0_mytype
    twothrd = 2.0_mytype/3.0_mytype
    thrfrth = 0.75_mytype
    im = xsize(1)
    jm = xsize(2)
    km = xsize(3)
    ! first substep k1 (internal and BC are asynchronous)
    call calcEuler_internal(drho,drhu,drhv,drhw,dret,rho,u,v,w,ien,pre,istep)
    call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time)
    ! calculation of the solution at time t
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm
        do i=1,im
          rhoOld(i,j,k) = rho(i,j,k)
          rhuOld(i,j,k) = rhoOld(i,j,k)*u(i,j,k)
          rhvOld(i,j,k) = rhoOld(i,j,k)*v(i,j,k)
          rhwOld(i,j,k) = rhoOld(i,j,k)*w(i,j,k)
          retOld(i,j,k) = rhoOld(i,j,k)*(ien(i,j,k) &
                        + 0.5_mytype*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2))
        enddo
      enddo
    enddo
    !$acc wait
    call calcEuler_atBoundary(drho,drhu,drhv,drhw,dret,rho,u,v,w,ien,pre,istep)
    call calcViscFlux(drho,drhu,drhv,drhw,dret,rho,u,v,w,ien,pre,tem,mu,ka)
    ! update of the solution with the first substep
    !$acc parallel loop collapse(3) default(present)
    do k=1,km
      do j=1,jm
        do i=1,im
          rho(i,j,k) =  rhoOld(i,j,k) + dt*drho(i,j,k)
          u(i,j,k)   = (rhuOld(i,j,k) + dt*drhu(i,j,k))/rho(i,j,k)
          v(i,j,k)   = (rhvOld(i,j,k) + dt*drhv(i,j,k))/rho(i,j,k)
          w(i,j,k)   = (rhwOld(i,j,k) + dt*drhw(i,j,k))/rho(i,j,k)
          ien(i,j,k) = (retOld(i,j,k) + dt*dret(i,j,k))/rho(i,j,k) &
                     - 0.5_mytype*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
        enddo
      enddo
    enddo
    ! update of the secondary variables
    call calcState_re(rho,ien,pre,tem,mu,ka,1,im,1,jm,1,km)

    ! second substep k2 (internal and BC are asynchronous)
    call calcEuler_internal(drho,drhu,drhv,drhw,dret,rho,u,v,w,ien,pre,istep)
    call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time+0.5*dt)
    !$acc wait
    call calcEuler_atBoundary(drho,drhu,drhv,drhw,dret,rho,u,v,w,ien,pre,istep)
    call calcViscFlux(drho,drhu,drhv,drhw,dret,rho,u,v,w,ien,pre,tem,mu,ka)
    ! update of the solution with first and second substeps
    !$acc parallel loop collapse(3) default(present) private(rho1,u1,v1,w1,ien1)
    do k=1,km
      do j=1,jm
        do i=1,im
          rho1 = rho(i,j,k)
          u1 = u(i,j,k)
          v1 = v(i,j,k)
          w1 = w(i,j,k)
          ien1 = ien(i,j,k)
          rho(i,j,k) =  thrfrth*rhoOld(i,j,k) + onefrth*(rho1    + dt*drho(i,j,k))
          u(i,j,k)   = (thrfrth*rhuOld(i,j,k) + onefrth*(rho1*u1 + dt*drhu(i,j,k)))/rho(i,j,k)
          v(i,j,k)   = (thrfrth*rhvOld(i,j,k) + onefrth*(rho1*v1 + dt*drhv(i,j,k)))/rho(i,j,k)
          w(i,j,k)   = (thrfrth*rhwOld(i,j,k) + onefrth*(rho1*w1 + dt*drhw(i,j,k)))/rho(i,j,k)
          ien(i,j,k) = (thrfrth*retOld(i,j,k) + onefrth*(rho1*(ien1 + 0.5_mytype*(u1**2+v1**2+w1**2)) &
                        + dt*dret(i,j,k)))/rho(i,j,k) - 0.5_mytype*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)                                           
        enddo
      enddo
    enddo
    ! update of the secondary variables
    call calcState_re(rho,ien,pre,tem,mu,ka,1,im,1,jm,1,km)

    ! final step (internal and BC are asynchronous)
    call calcEuler_internal(drho,drhu,drhv,drhw,dret,rho,u,v,w,ien,pre,istep)
    call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time+dt)
    !$acc wait
    call calcEuler_atBoundary(drho,drhu,drhv,drhw,dret,rho,u,v,w,ien,pre,istep)
    call calcViscFlux(drho,drhu,drhv,drhw,dret,rho,u,v,w,ien,pre,tem,mu,ka)
    ! final solution at time t+1
    !$acc parallel loop collapse(3) default(present) private(rho2,u2,v2,w2,ien2)
    do k=1,km
      do j=1,jm
        do i=1,im
          rho2 = rho(i,j,k)
          u2 = u(i,j,k)
          v2 = v(i,j,k)
          w2 = w(i,j,k)
          ien2 = ien(i,j,k)
          rho(i,j,k) =  onethrd*rhoOld(i,j,k) + twothrd*(rho2    + dt*drho(i,j,k))
          u(i,j,k)   = (onethrd*rhuOld(i,j,k) + twothrd*(rho2*u2 + dt*drhu(i,j,k)))/rho(i,j,k)
          v(i,j,k)   = (onethrd*rhvOld(i,j,k) + twothrd*(rho2*v2 + dt*drhv(i,j,k)))/rho(i,j,k)
          w(i,j,k)   = (onethrd*rhwOld(i,j,k) + twothrd*(rho2*w2 + dt*drhw(i,j,k)))/rho(i,j,k)
          ien(i,j,k) = (onethrd*retOld(i,j,k) + twothrd*(rho2*(ien2 + 0.5_mytype*(u2**2+v2**2+w2**2)) &
                        + dt*dret(i,j,k)))/rho(i,j,k) - 0.5_mytype*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
        enddo
      enddo
    enddo
    ! final update of the secondary variables
    call calcState_re(rho,ien,pre,tem,mu,ka,1,im,1,jm,1,km)
  end subroutine


end module mod_solve
