! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_solve
  use decomp_2d
  implicit none
  real(mytype), allocatable, dimension(:,:,:) :: drho1,drhu1,drhv1,drhw1,dret1, &
                                                 drho2,drhu2,drhv2,drhw2,dret2, &
                                                 drho3,drhu3,drhv3,drhw3,dret3, &
                                                 rhoOld,rhuOld,rhvOld,rhwOld,retOld
contains
! initialization Runge-Kutta
  subroutine init_rk3()
    use decomp_2d
    use mod_param
    implicit none
    ! 1st
    allocate( drho1(xsize(1), xsize(2), xsize(3)) )
    allocate( drhu1(xsize(1), xsize(2), xsize(3)) )
    allocate( drhv1(xsize(1), xsize(2), xsize(3)) )
    allocate( drhw1(xsize(1), xsize(2), xsize(3)) )
    allocate( dret1(xsize(1), xsize(2), xsize(3)) )
    ! 2nd
    allocate( drho2(xsize(1), xsize(2), xsize(3)) )
    allocate( drhu2(xsize(1), xsize(2), xsize(3)) )
    allocate( drhv2(xsize(1), xsize(2), xsize(3)) )
    allocate( drhw2(xsize(1), xsize(2), xsize(3)) )
    allocate( dret2(xsize(1), xsize(2), xsize(3)) )
    ! 3rd
    allocate( drho3(xsize(1), xsize(2), xsize(3)) )
    allocate( drhu3(xsize(1), xsize(2), xsize(3)) )
    allocate( drhv3(xsize(1), xsize(2), xsize(3)) )
    allocate( drhw3(xsize(1), xsize(2), xsize(3)) )
    allocate( dret3(xsize(1), xsize(2), xsize(3)) )
    ! old timestep
    allocate(rhoOld(xsize(1), xsize(2), xsize(3)) )
    allocate(rhuOld(xsize(1), xsize(2), xsize(3)) )
    allocate(rhvOld(xsize(1), xsize(2), xsize(3)) )
    allocate(rhwOld(xsize(1), xsize(2), xsize(3)) )
    allocate(retOld(xsize(1), xsize(2), xsize(3)) )
    !$acc enter data create(drho1,drhu1,drhv1,drhw1,dret1)
    !$acc enter data create(drho2,drhu2,drhv2,drhw2,dret2)
    !$acc enter data create(drho3,drhu3,drhv3,drhw3,dret3)
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
    do k=2,xsize(3)-1
      do j=1,xsize(2)
        do i=2,xsize(1)-1
          delx = (x(i+1)-x(i-1))/2.0
          delz = (z(k+1)-z(k-1))/2.0
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
      dt = CFL/max(dtConvInv, dtViscInv, dtCondInv)
      call mpi_allreduce(dt, tmp, 1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
      dt = tmp
    endif
    ! timestep with perturbation, a new CFL is calculated from dt_FFT (perturbation)
    if ((pert_calc==1) .OR. (BC_inl == "inlet_lst")) then
      CFL_new=dt_FFT * max(dtConvInv, dtViscInv, dtCondInv)
      call mpi_allreduce(CFL_new, tmp, 1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
      CFL_new=tmp
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
    real(mytype) :: dt,time
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl
    TYPE (DECOMP_INFO), intent(IN) :: part
    im = xsize(1)
    jm = xsize(2)
    km = xsize(3)
    ! first substep k1 (internal and BC are asynchronous)
    call calcEuler_internal(drho1,drhu1,drhv1,drhw1,dret1,rho,u,v,w,ien,pre,istep)
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
                        + 0.5*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2))
        enddo
      enddo
    enddo
    !$acc wait
    call calcEuler_atBoundary(drho1,drhu1,drhv1,drhw1,dret1,rho,u,v,w,ien,pre,istep)
    call calcViscFlux(drho1,drhu1,drhv1,drhw1,dret1,rho,u,v,w,ien,pre,tem,mu,ka)
    ! update of the solution with the first substep
    !$acc parallel loop collapse(3) default(present)
    do k=1,km
      do j=1,jm
        do i=1,im
          rho(i,j,k) =  rhoOld(i,j,k) + dt*0.5*drho1(i,j,k)
          u(i,j,k)   = (rhuOld(i,j,k) + dt*0.5*drhu1(i,j,k))/rho(i,j,k)
          v(i,j,k)   = (rhvOld(i,j,k) + dt*0.5*drhv1(i,j,k))/rho(i,j,k)
          w(i,j,k)   = (rhwOld(i,j,k) + dt*0.5*drhw1(i,j,k))/rho(i,j,k)
          ien(i,j,k) = (retOld(i,j,k) + dt*0.5*dret1(i,j,k))/rho(i,j,k) &
                     - 0.5*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
        enddo
      enddo
    enddo
    ! update of the secondary variables
    call calcState_re(rho,ien,pre,tem,mu,ka,1,im,1,jm,1,km)
    ! second substep k2 (internal and BC are asynchronous)
    call calcEuler_internal(drho2,drhu2,drhv2,drhw2,dret2,rho,u,v,w,ien,pre,istep)
    call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time+0.5*dt)
    !$acc wait
    call calcEuler_atBoundary(drho2,drhu2,drhv2,drhw2,dret2,rho,u,v,w,ien,pre,istep)
    call calcViscFlux(drho2,drhu2,drhv2,drhw2,dret2,rho,u,v,w,ien,pre,tem,mu,ka)
    ! update of the solution with first and second substeps
    !$acc parallel loop collapse(3) default(present)
    do k=1,km
      do j=1,jm
        do i=1,im
          rho(i,j,k) =  rhoOld(i,j,k) + dt*(-drho1(i,j,k)+2.0*drho2(i,j,k))
          u(i,j,k)   = (rhuOld(i,j,k) + dt*(-drhu1(i,j,k)+2.0*drhu2(i,j,k)))/rho(i,j,k)
          v(i,j,k)   = (rhvOld(i,j,k) + dt*(-drhv1(i,j,k)+2.0*drhv2(i,j,k)))/rho(i,j,k)
          w(i,j,k)   = (rhwOld(i,j,k) + dt*(-drhw1(i,j,k)+2.0*drhw2(i,j,k)))/rho(i,j,k)
          ien(i,j,k) = (retOld(i,j,k) + dt*(-dret1(i,j,k)+2.0*dret2(i,j,k)))/rho(i,j,k) & 
                     - 0.5*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)                                           
        enddo
      enddo
    enddo
    ! update of the secondary variables
    call calcState_re(rho,ien,pre,tem,mu,ka,1,im,1,jm,1,km)
    ! second substep k2 (internal and BC are asynchronous)
    call calcEuler_internal(drho3,drhu3,drhv3,drhw3,dret3,rho,u,v,w,ien,pre,istep)
    call setBC(part,rho,u,v,w,ien,pre,tem,mu,ka,rho_bl,u_bl,v_bl,w_bl,ien_bl,pre_bl,tem_bl,mu_bl,ka_bl,time+dt)
    !$acc wait
    call calcEuler_atBoundary(drho3,drhu3,drhv3,drhw3,dret3,rho,u,v,w,ien,pre,istep)
    call calcViscFlux(drho3,drhu3,drhv3,drhw3,dret3,rho,u,v,w,ien,pre,tem,mu,ka)
    ! final solution at time t+1
    !$acc parallel loop collapse(3) default(present)
    do k=1,km
      do j=1,jm
        do i=1,im
          rho(i,j,k) =  rhoOld(i,j,k) + dt*(drho1(i,j,k)+4.0*drho2(i,j,k)+drho3(i,j,k))/6.0
          u(i,j,k)   = (rhuOld(i,j,k) + dt*(drhu1(i,j,k)+4.0*drhu2(i,j,k)+drhu3(i,j,k))/6.0)/rho(i,j,k)
          v(i,j,k)   = (rhvOld(i,j,k) + dt*(drhv1(i,j,k)+4.0*drhv2(i,j,k)+drhv3(i,j,k))/6.0)/rho(i,j,k)
          w(i,j,k)   = (rhwOld(i,j,k) + dt*(drhw1(i,j,k)+4.0*drhw2(i,j,k)+drhw3(i,j,k))/6.0)/rho(i,j,k)
          ien(i,j,k) = (retOld(i,j,k) + dt*(dret1(i,j,k)+4.0*dret2(i,j,k)+dret3(i,j,k))/6.0)/rho(i,j,k) &
                     - 0.5*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
        enddo
      enddo
    enddo
    ! final update of the secondary variables
    call calcState_re(rho,ien,pre,tem,mu,ka,1,im,1,jm,1,km)
  end subroutine
end module mod_solve
