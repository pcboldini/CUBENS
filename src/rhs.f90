! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! Right-Hand-Side module

module mod_rhs
  use decomp_2d
  use mod_param
  use mod_grid
  use mod_finitediff
  use mod_timer
  implicit none
  ! define allocation RHS tensors, dilatation, and sponge parameters
  real(mytype), allocatable, dimension(:,:,:)   :: dil
  real(mytype), allocatable, dimension(:,:,:)   :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e
  real(mytype), allocatable, dimension(:,:)     :: r_ref, ru_ref, rv_ref, rw_ref, ret_ref, p_ref
  real(mytype), allocatable, dimension(:)       :: spSigX, spSigZ
  real(mytype), allocatable, dimension(:,:,:)   :: rhoe
contains


! initialization allocated variables and sponge
  subroutine init_rhs()
    use decomp_2d
    use mod_param
    use mod_grid
    use mod_init
    implicit none
    integer :: i,j,k,kk
    real(mytype) :: fz
    ! allocation
    allocate(          dil(1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate( rhs_r(xsize(1), xsize(2), xsize(3)) )
    allocate( rhs_u(xsize(1), xsize(2), xsize(3)) )
    allocate( rhs_v(xsize(1), xsize(2), xsize(3)) )
    allocate( rhs_w(xsize(1), xsize(2), xsize(3)) )
    allocate( rhs_e(xsize(1), xsize(2), xsize(3)) )
    allocate( rhoe(1-nHalo:xsize(1)+nHalo, 1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo) )
    allocate(  p_ref(xsize(1),xsize(3)))
    ! allocate sponge if needed
    if ((spInlLen.gt.0.0).or.(spOutLen.gt.0.0).or.(spTopLen.gt.0.0)) then 
      allocate(  r_ref(xsize(1),xsize(3)))
      allocate( ru_ref(xsize(1),xsize(3)))
      allocate( rv_ref(xsize(1),xsize(3)))
      allocate( rw_ref(xsize(1),xsize(3)))
      allocate(ret_ref(xsize(1),xsize(3)))
      allocate(spSigX(xsize(1)))
      allocate(spSigZ(xsize(3)))
      ! set sponge damping coefficient in X direction
      spSigX = 0.0_mytype
      do i=1, xsize(1)
        if ((spTopLen .gt. 0.0_mytype) .and. (x(i) .ge. (len_x-spTopLen))) then 
          spSigX(i) = spTopStr*((x(i) - (len_x-spTopLen))/spTopLen)**spTopExp
        endif
      enddo
      ! set sponge damping coefficient in Z direction
      spSigZ = 0.0_mytype
      spInLen_Reend=(Re)**0.5*( (spInlLen+zStartDNS)**0.5 )
      spOutLen_Resta=(Re)**0.5*( (zEndDNS-spOutLen)**0.5 )
      do k=1, xsize(3)
        if ((spInlLen .gt. 0.0) .and. (z(k) .le. spInlLen)) then
          spSigZ(k) = spInlStr*((spInlLen - z(k))/spInlLen)**spInlExp
        endif
        if ((spOutLen .gt. 0.0) .and. (z(k) .ge. (len_z-spOutLen))) then
          spSigZ(k) = spOutStr*((z(k) - (len_z-spOutLen))/spOutLen)**spOutExp
        endif
      enddo
    endif
    !$acc enter data copyin(rhoe,dil) async 
    !$acc enter data copyin(spSigX,spSigZ) async
  end subroutine


! calculation of divergence
  subroutine calcDivergence(u,v,w)
    use decomp_2d
    implicit none
    integer i,j,k,c
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: u,v,w
    !$acc parallel loop collapse(3) default(present)
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          dil(i,j,k) = 0.0_mytype
          !$acc loop seq
          do c = 1, nStencilVisc
            dil(i,j,k) = dil(i,j,k) &
                        + visc_ddx(c)*(u(i+c,j,k) - u(i-c,j,k))*xp(i)  &
                        + visc_ddy(c)*(v(i,j+c,k) - v(i,j-c,k))        &
                        + visc_ddz(c)*(w(i,j,k+c) - w(i,j,k-c))*zp(k)
          enddo
        enddo
      enddo
    enddo
  end subroutine


! calculation of convective fluxes at internal mesh points
  subroutine calcEuler_internal(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e,p,istep) 
    use decomp_2d
    use mod_param
    use mod_halo
    use mod_boundary
    implicit none
    integer :: im,jm,km, i1,i2,k1,k2,istep
    integer :: i,j,k
    real(mytype), dimension(:,:,:) :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:)  :: r,u,v,w,e,p
    real(mytype) :: dprex, dprey, dprez
    integer :: s
    real(mytype) :: rhoip, rhoim, rhojp, rhojm, rhokp, rhokm
    real(mytype) :: uip, uim, ujp, ujm, ukp, ukm
    real(mytype) :: vip, vim, vjp, vjm, vkp, vkm
    real(mytype) :: wip, wim, wjp, wjm, wkp, wkm
    real(mytype) :: eip, eim, ejp, ejm, ekp, ekm
    real(mytype) :: rhoeip, rhoeim, rhoejp, rhoejm, rhoekp, rhoekm
    real(mytype) :: rhoa, ua, va, wa, pa, ea, rhoea, za, xa

    im = xsize(1)
    jm = xsize(2)
    km = xsize(3)
    ! initializiation with null
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,km
      do j=1,jm
        do i=1,im
          rhs_r(i,j,k) = 0.0_mytype
          rhs_u(i,j,k) = 0.0_mytype
          rhs_v(i,j,k) = 0.0_mytype
          rhs_w(i,j,k) = 0.0_mytype
          rhs_e(i,j,k) = 0.0_mytype
        enddo
      enddo
    enddo
    if ( keep_flag == "pep" ) then
      ! calculation of rho*e, needed for internal energy flux
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,km
        do j=1,jm
          do i=1,im
            rhoe(i,j,k)=r(i,j,k)*e(i,j,k)
          enddo
        enddo
      enddo
    endif
    ! only in the internal cells
    if ((im > 2*nHalo) .and. (jm > 2*nHalo) .and. (km > 2*nHalo)) then
      !$acc parallel loop collapse(3) default(present) private(rhoa,ua,va,wa,pa,rhoea,za,xa) async(1)
      do k=1+nHalo,km-nHalo 
        do j=1+nHalo,jm-nHalo
          do i=1+nHalo,im-nHalo
            ! routine for fluxes calculation 
            if ( keep_flag == "classic" ) then 
#include "eulerFlux.f90"
            elseif ( keep_flag == "pep" ) then 
#include "eulerFlux_pep.f90"
            endif
          enddo
        enddo
      enddo
    endif
  end subroutine


! calculation of convective fluxes at boundary mesh points, top and bottom boundaries excluded
  subroutine calcEuler_atBoundary(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e,p,istep)
    use decomp_2d
    use mod_param
    use mod_halo
    use mod_boundary
    implicit none
    integer :: im,jm,km,i1,i2,k1,k2,istep,i,j,k,iminHalo,imaxHalo,jminHalo,jmaxHalo,kminHalo,kmaxHalo
    real(mytype), dimension(:,:,:) :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:)  :: r,u,v,w,e,p
    real(mytype) :: dprex, dprey, dprez
    integer :: s
    real(mytype) :: rhoip, rhoim, rhojp, rhojm, rhokp, rhokm
    real(mytype) :: uip, uim, ujp, ujm, ukp, ukm
    real(mytype) :: vip, vim, vjp, vjm, vkp, vkm
    real(mytype) :: wip, wim, wjp, wjm, wkp, wkm
    real(mytype) :: eip, eim, ejp, ejm, ekp, ekm
    real(mytype) :: rhoeip, rhoeim, rhoejp, rhoejm, rhoekp, rhoekm
    real(mytype) :: rhoa, ua, va, wa, pa, ea, rhoea, za, xa

    im = xsize(1)
    jm = xsize(2)
    km = xsize(3)
    ! boundary cells 
    i1 = 1
    i2 = im 
    if (perBC(1) .eqv. .false.) then 
      i1 = 2
      i2 = im-1 
    endif 
    k1 = 1 
    k2 = km
    if ((perBC(3) .eqv. .false.) .and. (neigh%inlet )) then 
      k1 = 2
    endif
    if ((perBC(3) .eqv. .false.) .and. (neigh%outlet)) then
      k2 = km-1
    endif
    imaxHalo=min(nHalo,i2)
    iminHalo=max(im-nHalo+1,nHalo+1)
    jmaxHalo=min(nHalo,jm)
    jminHalo=max(jm-nHalo+1,nHalo+1)
    kmaxHalo=min(nHalo,k2)
    kminHalo=max(km-nHalo+1,nHalo+1)
    if ( keep_flag == "pep" ) then
      ! calculation of rho*e, needed for internal energy flux
      ! k face min
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1-nHalo,0
        do j=1,jm
          do i=1,im
            rhoe(i,j,k)=r(i,j,k)*e(i,j,k)
          enddo
        enddo
      enddo
      ! k face max
      !$acc parallel loop collapse(3) default(present) async(2)
      do k=km+1,km+nHalo
        do j=1,jm
          do i=1,im
            rhoe(i,j,k)=r(i,j,k)*e(i,j,k)
          enddo
        enddo
      enddo
      ! j face min
      !$acc parallel loop collapse(3) default(present) async(3)
      do k=1,km
        do j=1-nHalo,0
          do i=1,im
            rhoe(i,j,k)=r(i,j,k)*e(i,j,k)
          enddo
        enddo
      enddo
      ! j face max
      !$acc parallel loop collapse(3) default(present) async(4)
      do k=1,km
        do j=jm+1,jm+nHalo
          do i=1,im
            rhoe(i,j,k)=r(i,j,k)*e(i,j,k)
          enddo
        enddo
      enddo
      ! i face min
      !$acc parallel loop collapse(3) default(present) async(5)
      do k=1,km
        do j=1,jm
          do i=1-nHalo,0
            rhoe(i,j,k)=r(i,j,k)*e(i,j,k)
          enddo
        enddo
      enddo
      ! i face max
      !$acc parallel loop collapse(3) default(present) async(6)
      do k=1,km
        do j=1,jm
          do i=im+1,im+nHalo
            rhoe(i,j,k)=r(i,j,k)*e(i,j,k)
          enddo
        enddo
      enddo
      !$acc wait
    endif
    ! loop over k-minus boundary
    !$acc parallel loop collapse(3) default(present) private(rhoa,ua,va,wa,pa,ea,rhoea,za,xa) async(1)
    do k=k1,kmaxHalo
      do j=1,jm
        do i=i1,i2
          if ( keep_flag == "classic" ) then 
#include "eulerFlux.f90"
          elseif ( keep_flag == "pep" ) then 
#include "eulerFlux_pep.f90"
          endif
        enddo
      enddo
    enddo
    ! loop over k-plus boundary
    if (km>nHalo) then
      !$acc parallel loop collapse(3) default(present) private(rhoa,ua,va,wa,pa,ea,rhoea,za,xa) async(2)
      do k=kminHalo,k2
        do j=1,jm
          do i=i1,i2
            if ( keep_flag == "classic" ) then 
#include "eulerFlux.f90"
            elseif ( keep_flag == "pep" ) then 
#include "eulerFlux_pep.f90"
            endif
          enddo
        enddo
      enddo
    endif
    ! loop over j-minus boundary 
    !$acc parallel loop collapse(3) default(present) private(rhoa,ua,va,wa,pa,ea,rhoea,za,xa) async(3)
    do k=kmaxHalo+1,kminHalo-1
      do j=1,jmaxHalo
        do i=i1,i2
          if ( keep_flag == "classic" ) then 
#include "eulerFlux.f90"
          elseif ( keep_flag == "pep" ) then 
#include "eulerFlux_pep.f90"
          endif
        enddo
      enddo
    enddo
    ! loop over j-plus boundary 
    if (jm>nHalo) then
      !$acc parallel loop collapse(3) default(present) private(rhoa,ua,va,wa,pa,ea,rhoea,za,xa) async(4)
      do k=kmaxHalo+1,kminHalo-1
        do j=jminHalo,jm
          do i=i1,i2
            if ( keep_flag == "classic" ) then 
#include "eulerFlux.f90"
            elseif ( keep_flag == "pep" ) then 
#include "eulerFlux_pep.f90"
            endif
          enddo
        enddo
      enddo
    endif   
    if (jm > 2*nHalo) then
      ! loop over i-minus boundary 
      !$acc parallel loop collapse(3) default(present) private(rhoa,ua,va,wa,pa,ea,rhoea,za,xa) async(5)
      do k=kmaxHalo+1,kminHalo-1
        do j=jmaxHalo+1,jminHalo-1
          do i=i1,imaxHalo
            if ( keep_flag == "classic" ) then 
#include "eulerFlux.f90"
            elseif ( keep_flag == "pep" ) then 
#include "eulerFlux_pep.f90"
            endif
          enddo
        enddo
      enddo
      ! loop over i-plus boundary 
      if (im>nHalo) then
        !$acc parallel loop collapse(3) default(present) private(rhoa,ua,va,wa,pa,ea,rhoea,za,xa) async(6)
        do k=kmaxHalo+1,kminHalo-1
          do j=jmaxHalo+1,jminHalo-1
            do i=iminHalo,i2
              if ( keep_flag == "classic" ) then 
#include "eulerFlux.f90"
              elseif ( keep_flag == "pep" ) then 
#include "eulerFlux_pep.f90"
              endif
            enddo
          enddo
        enddo
      endif 
    endif
    !$acc wait
  end subroutine


! call of the viscous (diffusive) RHS at all cells 
! call of the boundary conditions and sponge
  subroutine calcViscFlux(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e,p,tem,mu,ka)
    use decomp_2d
    use mod_param
    use mod_halo
    use mod_boundary
    implicit none
    real(mytype), dimension(:,:,:) :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:)  :: r,u,v,w,e,p,tem,mu,ka
    ! calculation of the viscous (diffusive) fluxes
    call calcRHSVisc(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e,tem,mu,ka)
    ! add streamwise pressure term
    !$acc kernels default(present)
    if (dpdz /= 0.0) then 
      rhs_w = rhs_w + dpdz
      rhs_e = rhs_e + dpdz*w(1:xsize(1),1:xsize(2),1:xsize(3))
    endif
    !$acc end kernels
    ! bottom boundary conditions
    if (perBC(1) .eqv. .false.) then
      call setBC_RHS_Bot(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e,p)
    endif
    ! top boundary conditions
    if (perBC(1) .eqv. .false.) then
      call setBC_RHS_Top(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e,p)
    endif
    ! inlet boundary conditions
    if ((perBC(3) .eqv. .false.) .and. (neigh%inlet)) then
      call setBC_RHS_Inl(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e,p)
    endif
    ! outlet boundary conditions
    if ((perBC(3) .eqv. .false.) .and. (neigh%outlet)) then
      call setBC_RHS_Out(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e,p,p_ref)
    endif
    !$acc wait    
    ! sponge enforcement 
    if ((spInlLen.gt.0.0).or.(spOutLen.gt.0.0).or.(spTopLen.gt.0.0)) then 
      call sponge(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e)
    endif
  end subroutine


! calculation of the viscous (diffusive) RHS at all cells 
  subroutine calcRHSVisc(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e, r,u,v,w,e,T,mu,ka) 
    use decomp_2d
    use mod_param
    use mod_halo
    use mod_boundary
    use mod_finitediff
    implicit none
    integer  :: i,j,k,c
    real(mytype), dimension(:,:,:), intent(inout   )  :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:), intent(inout   )  :: r,u,v,w,e,T,mu,ka
    real(mytype) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, &
                    d2ux, d2vx, d2wx, d2uy, d2vy, d2wy, d2uz, d2vz, d2wz, &
                    dTx, dTy, dTz, d2Tx, d2Ty, d2Tz, &
                    ddilx, ddily, ddilz, &
                    sxx, sxy, sxz, syy, syz, szz, &
                    dtauxj, dtauyj, dtauzj, & 
                    dmux, dmuy, dmuz, dkax, dkay, dkaz
    real(mytype) :: ua, va, wa, Ta, dila, mua, kaa, dilid, diljd, dilkd
    real(mytype) :: uip, uim, ujp, ujm, ukp, ukm
    real(mytype) :: vip, vim, vjp, vjm, vkp, vkm
    real(mytype) :: wip, wim, wjp, wjm, wkp, wkm
    real(mytype) :: Tip, Tim, Tjp, Tjm, Tkp, Tkm
    ! calculation of the divergence
    call calcDivergence(u,v,w)
    !$acc parallel loop collapse(3) default(present)
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          dil(i,j,k) = dil(i,j,k) / 3
        enddo
      enddo
    enddo
    ! MPI update for the divergence
#if defined(_OPENACC)
    call haloUpdateMult_i((/.true.,.true.,.true./), xsize, dil) 
    call haloUpdate1_jk((/.true.,.true.,.true./), xsize, dil)
#else
    call haloUpdateMult_CPU((/.true.,.true.,.true./), xsize, dil)  
#endif
    ! apply boundary conditions on divergence
    if (perBC(1) .eqv. .false.) then
      !$acc parallel loop collapse(2) default(present)
      do k=1,xsize(3)
        do j=1,xsize(2)
          !$acc loop seq
          do i=-1,0
            dil(i,j,k) = dil(2-i,j,k)
          enddo
          !$acc loop seq
          do i=1,2
            dil(xsize(1)+i,j,k) = dil(xsize(1)-i,j,k)
          enddo
        enddo
      enddo
    endif
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dtauzj,dtauyj,dtauxj,dmuz,szz,syz,syy,sxz,sxy,sxx) &
    !$acc private(d2Tz,d2Ty,d2Tx,d2wz,d2wy,d2wx,d2vz,d2vy,d2vx,d2uz,d2uy,d2ux) &
    !$acc private(dkaz,dkay,dkax,ddilz,dmuz,dmuy,dmux,dTz,dTy,dTx,dwz,dwy,dwx,dvz,dvy,dvx,duz,duy,dux,ddilz,ddily,ddilx) &
    !$acc private(ua,va,wa,Ta,dila,mua,kaa,dilid,diljd,dilkd) &
    !$acc private(uip,uim,ujp,ujm,ukp,ukm,vip,vim,vjp,vjm,vkp,vkm,wip,wim,wjp,wjm,wkp,wkm,Tip,Tim,Tjp,Tjm,Tkp,Tkm)
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          ua    = u(i,j,k)
          va    = v(i,j,k)
          wa    = w(i,j,k)
          Ta    = T(i,j,k)
          dila  = dil(i,j,k)
          mua   = mu(i,j,k)
          kaa   = ka(i,j,k)
          ddilx = 0.0_mytype
          ddily = 0.0_mytype
          ddilz = 0.0_mytype
          dux   = 0.0_mytype
          duy   = 0.0_mytype
          duz   = 0.0_mytype
          dvx   = 0.0_mytype
          dvy   = 0.0_mytype
          dvz   = 0.0_mytype
          dwx   = 0.0_mytype
          dwy   = 0.0_mytype
          dwz   = 0.0_mytype
          dTx   = 0.0_mytype
          dTy   = 0.0_mytype
          dTz   = 0.0_mytype
          dmux  = 0.0_mytype
          dmuy  = 0.0_mytype
          dmuz  = 0.0_mytype
          dkax  = 0.0_mytype
          dkay  = 0.0_mytype
          dkaz  = 0.0_mytype
          d2ux = visc_d2dx2(0,i)*ua
          d2uy = visc_d2dy2(0)  *ua
          d2uz = visc_d2dz2(0,k)*ua
          d2vx = visc_d2dx2(0,i)*va
          d2vy = visc_d2dy2(0)  *va
          d2vz = visc_d2dz2(0,k)*va
          d2wx = visc_d2dx2(0,i)*wa
          d2wy = visc_d2dy2(0)  *wa
          d2wz = visc_d2dz2(0,k)*wa
          d2Tx = visc_d2dx2(0,i)*Ta
          d2Ty = visc_d2dy2(0)  *Ta
          d2Tz = visc_d2dz2(0,k)*Ta
          !$acc loop seq
          do c = 1, nStencilVisc
            ! u
            uip   = u(i+c,j,k)
            uim   = u(i-c,j,k)
            ujp   = u(i,j+c,k)
            ujm   = u(i,j-c,k)
            ukp   = u(i,j,k+c)
            ukm   = u(i,j,k-c)
            ! v
            vip   = v(i+c,j,k)
            vim   = v(i-c,j,k)
            vjp   = v(i,j+c,k)
            vjm   = v(i,j-c,k)
            vkp   = v(i,j,k+c)
            vkm   = v(i,j,k-c)
            ! w
            wip   = w(i+c,j,k)
            wim   = w(i-c,j,k)
            wjp   = w(i,j+c,k)
            wjm   = w(i,j-c,k)
            wkp   = w(i,j,k+c)
            wkm   = w(i,j,k-c)
            ! T
            Tip   = T(i+c,j,k)
            Tim   = T(i-c,j,k)
            Tjp   = T(i,j+c,k)
            Tjm   = T(i,j-c,k)
            Tkp   = T(i,j,k+c)
            Tkm   = T(i,j,k-c)
            ! divergence
            dilid = dil(i+c,j,k) - dil(i-c,j,k)
            diljd = dil(i,j+c,k) - dil(i,j-c,k)
            dilkd = dil(i,j,k+c) - dil(i,j,k-c)
            ! first-order derivatives
            ddilx = ddilx + visc_ddx(c)*(dilid)*xp(i)
            ddily = ddily + visc_ddy(c)*(diljd)
            ddilz = ddilz + visc_ddz(c)*(dilkd)*zp(k)
            dux   = dux   + visc_ddx(c)*( uip - uim )*xp(i)
            duy   = duy   + visc_ddy(c)*( ujp - ujm )
            duz   = duz   + visc_ddz(c)*( ukp - ukm )*zp(k)
            dvx   = dvx   + visc_ddx(c)*( vip - vim )*xp(i)
            dvy   = dvy   + visc_ddy(c)*( vjp - vjm )
            dvz   = dvz   + visc_ddz(c)*( vkp - vkm )*zp(k)
            dwx   = dwx   + visc_ddx(c)*( wip - wim )*xp(i)
            dwy   = dwy   + visc_ddy(c)*( wjp - wjm )
            dwz   = dwz   + visc_ddz(c)*( wkp - wkm )*zp(k)
            dTx   = dTx   + visc_ddx(c)*( Tip - Tim )*xp(i)
            dTy   = dTy   + visc_ddy(c)*( Tjp - Tjm )
            dTz   = dTz   + visc_ddz(c)*( Tkp - Tkm )*zp(k)
            dmux  = dmux  + visc_ddx(c)*(mu(i+c,j,k) - mu(i-c,j,k))*xp(i)
            dmuy  = dmuy  + visc_ddy(c)*(mu(i,j+c,k) - mu(i,j-c,k))
            dmuz  = dmuz  + visc_ddz(c)*(mu(i,j,k+c) - mu(i,j,k-c))*zp(k)
            dkax  = dkax  + visc_ddx(c)*(ka(i+c,j,k) - ka(i-c,j,k))*xp(i)
            dkay  = dkay  + visc_ddy(c)*(ka(i,j+c,k) - ka(i,j-c,k))
            dkaz  = dkaz  + visc_ddz(c)*(ka(i,j,k+c) - ka(i,j,k-c))*zp(k)
            ! second-order derivatives
            d2ux = d2ux + visc_d2dx2(c,i)*  uip + visc_d2dx2(-c,i)*uim 
            d2uy = d2uy + visc_d2dy2(c)  * (ujp +                  ujm )
            d2uz = d2uz + visc_d2dz2(c,k)*  ukp + visc_d2dz2(-c,k)*ukm  
            d2vx = d2vx + visc_d2dx2(c,i)*  vip + visc_d2dx2(-c,i)*vim 
            d2vy = d2vy + visc_d2dy2(c)  * (vjp +                  vjm )
            d2vz = d2vz + visc_d2dz2(c,k)*  vkp + visc_d2dz2(-c,k)*vkm  
            d2wx = d2wx + visc_d2dx2(c,i)*  wip + visc_d2dx2(-c,i)*wim 
            d2wy = d2wy + visc_d2dy2(c)  * (wjp +                  wjm)
            d2wz = d2wz + visc_d2dz2(c,k)*  wkp + visc_d2dz2(-c,k)*wkm  
            d2Tx = d2Tx + visc_d2dx2(c,i)*  Tip + visc_d2dx2(-c,i)*Tim  
            d2Ty = d2Ty + visc_d2dy2(c)  * (Tjp +                  Tjm )
            d2Tz = d2Tz + visc_d2dz2(c,k)*  Tkp + visc_d2dz2(-c,k)*Tkm 
          enddo
          ! stress tensor
          sxx = 2.0*(dux - dila)
          sxy =      duy + dvx
          sxz =      duz + dwx
          syy = 2.0*(dvy - dila)
          syz =      dvz + dwy
          szz = 2.0*(dwz - dila)
          dtauxj = mua*(d2ux + d2uy + d2uz + ddilx) + dmux*sxx + dmuy*sxy + dmuz*sxz
          dtauyj = mua*(d2vx + d2vy + d2vz + ddily) + dmux*sxy + dmuy*syy + dmuz*syz
          dtauzj = mua*(d2wx + d2wy + d2wz + ddilz) + dmux*sxz + dmuy*syz + dmuz*szz
          ! RHS
          rhs_u(i,j,k) = rhs_u(i,j,k) + dtauxj
          rhs_v(i,j,k) = rhs_v(i,j,k) + dtauyj
          rhs_w(i,j,k) = rhs_w(i,j,k) + dtauzj
          rhs_e(i,j,k) = rhs_e(i,j,k) &
                       + ua*dtauxj + va*dtauyj + wa*dtauzj & 
                       + mua*(sxx*dux + sxy*duy + sxz*duz  &
                             +sxy*dvx + syy*dvy + syz*dvz  &
                             +sxz*dwx + syz*dwy + szz*dwz) & 
                       + ka(i,j,k)*(d2Tx + d2Ty + d2Tz) + dkax*dTx + dkay*dTy + dkaz*dTz
        enddo
      enddo
    enddo
  end subroutine
! calculation of the sponge
  subroutine sponge(rhs_r,rhs_u,rhs_v,rhs_w,rhs_e,r,u,v,w,e)
    use decomp_2d
    use mod_param
    implicit none
    integer :: i,j,k, im, jm, km
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: r,u,v,w,e
    real(mytype), dimension(:,:,:) :: rhs_r,rhs_u,rhs_v,rhs_w,rhs_e
    real(mytype)                   :: rhoa, spSigZa
    !$acc parallel loop collapse(3) default(present) private(rhoa,spSigZa)
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1) 
          spSigZa = spSigZ(k)
          rhoa        = r(i,j,k)
          rhs_r(i,j,k)=rhs_r(i,j,k)+(spSigX(i)+spSigZa)*(  r_ref(i,k)-rhoa         )
          rhs_u(i,j,k)=rhs_u(i,j,k)+(spSigX(i)+spSigZa)*( ru_ref(i,k)-rhoa*u(i,j,k))
          rhs_v(i,j,k)=rhs_v(i,j,k)+(spSigX(i)+spSigZa)*( rv_ref(i,k)-rhoa*v(i,j,k))
          rhs_w(i,j,k)=rhs_w(i,j,k)+(spSigX(i)+spSigZa)*( rw_ref(i,k)-rhoa*w(i,j,k))
          rhs_e(i,j,k)=rhs_e(i,j,k)+(spSigX(i)+spSigZa)*(ret_ref(i,k) & 
                                                       - rhoa*(e(i,j,k)+0.5*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) ))
        enddo
      enddo
    enddo
  end subroutine


! print initial parameters for RHS (sponge)
  subroutine print_init_rhs()
    use decomp_2d
    use mod_param
    use mod_grid
    use mod_finitediff
    implicit none
    real(mytype) :: ReBlasiusEnd
    ! calculation of the outlet Reynolds number
    ReBlasiusEnd=(Re)**0.5*( (zEndDNS)**0.5 )  
    if (nrank == 0) then
#if defined(BL)
      write(stdout,* ) 'Sponge'
      write(stdout,'(A)') 'o--------------------------------------------------o'
#endif
      if (spInlLen>0.0) then
        write(stdout,* ) 'Inlet:                             '
        write(stdout,"(A26,F10.4)") 'spInlLen from 0 to ',spInlLen
        write(stdout,"(A24,F10.4,A3,F10.4)") 'or Re_delta from ',Redelta_start,' to ', spInLen_Reend  
        write(stdout,"(A32,F10.4)") 'Max. damping coefficient=',spInlStr
        write(stdout,* )  
      else
      endif 
      if (spOutLen>0.0) then
        write(stdout,* ) 'Outlet:                             '
        write(stdout,"(A22,F10.4,A3,F10.4)") 'spOutLen from ',(len_z-spOutLen),' to ', len_z
        write(stdout,"(A25,F10.4,A3,F10.4)") 'or Re_delta from ',spOutLen_Resta,' to ', ReBlasiusEnd          
        write(stdout,"(A33,F10.4)") 'Max. damping coefficient=',spOutStr   
        write(stdout,* ) 
      else
      endif  
      if (spTopLen>0.0) then 
        write(stdout,* ) 'Free stream:                         '
        write(stdout,"(A22,F10.4,A3,F10.4)") 'spTopLen from ',(len_x-spTopLen),' to ', len_x
        write(stdout,"(A33,F10.4)") 'Max. damping coefficient=',spTopStr 
        write(stdout,'(A)') 'o--------------------------------------------------o'
        write(stdout,* )
      else  
      endif               
    endif
  end subroutine

  
end module
