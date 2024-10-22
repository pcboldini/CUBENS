! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! auxiliary routines for interpolation

module mod_interpolate
  use decomp_2d
  use mod_halo
  use mod_param
  use mod_auxl
  use mod_grid
  implicit none
  real(mytype), allocatable, dimension(:,:,:) :: rbo,ubo,vbo,wbo,ebo
  real(mytype), allocatable, dimension(:,:,:) :: rbn,ubn,vbn,wbn,ebn
  real(mytype), allocatable, dimension(:) :: xold, yold, zold, zold_global
  real(mytype), allocatable, dimension(:) :: xnew, ynew, znew, znew_global
  TYPE (DECOMP_INFO) :: part1, part2
  contains


! initializing the domain
  subroutine initDomain(xmesh_type,part,im,jm,km,Lx,Ly,Lz,&
                        ReTau,stretchx, &
                        zmesh_type,z_1,z_2,bumpz1,bumpz2,zplus_min,zplus_max, &
                        zpert_1, zpert_2, bumpzpert1, bumpzpert2, zpluspert_min, &
                        x,y,z,z_global,r,u,v,w,e)
    use decomp_2d
    implicit none 
    integer :: j,k,im,jm,km,ierr
    real(mytype) :: Lx, Ly, Lz, ReTau, stretchx, dy, z_1, z_2, bumpz1, bumpz2, zplus_min, zplus_max
    real(mytype) :: zpert_1, zpert_2, bumpzpert1, bumpzpert2, zpluspert_min
    real(mytype), allocatable, dimension(:)     :: x,y,z,z_global
    real(mytype), allocatable, dimension(:)     :: dummy1,dummy2,tmp
    real(mytype), allocatable, dimension(:,:,:) :: r,u,v,w,e
    character(len=30) :: xmesh_type, zmesh_type
    TYPE (DECOMP_INFO), intent(IN) :: part
    allocate( dummy1(im), dummy2(km), tmp(km))
    allocate( x(part%xsz(1)) )
    allocate( y(part%xsz(2)) )
    allocate( z(part%xsz(3)), z_global(km))
    ! x-mesh
#if !defined(CHA)
    ! for boundary layer and TGV  
    call xDistribution(xmesh_type, ReTau, stretchx, Lx, part%xsz(1), x, dummy1, dummy1)
#elif defined(CHA)
    ! for channel
    call xDistribution_CHA(stretchx, Lx, part%xsz(1), x, dummy1, dummy1)
#endif 
    ! y-mesh
    dy = Ly/(jm)
    call yDistribution(jm, y, dy)
    ! z-mesh
    call zDistribution(zmesh_type,part%xst(3),part%xsz(3), km, Lz, &
                        z_1, z_2, bumpz1, bumpz2, zplus_min, zplus_max, &
                        zpert_1, zpert_2, bumpzpert1, bumpzpert2, zpluspert_min, &
                        z, dummy2, dummy2, z_global, dummy2, dummy2)
    ! allocate primary variables
    allocate(r(part%xsz(1), part%xsz(2), part%xsz(3)))
    allocate(u(part%xsz(1), part%xsz(2), part%xsz(3)))
    allocate(v(part%xsz(1), part%xsz(2), part%xsz(3)))
    allocate(w(part%xsz(1), part%xsz(2), part%xsz(3)))
    allocate(e(part%xsz(1), part%xsz(2), part%xsz(3)))
  end subroutine


! initializing solution
  subroutine initSolution(part,x,z,r1,u1,v1,w1,e1)
    use mod_init
    real(mytype), dimension(:) :: x,z
    real(mytype), dimension(:,:,:) :: r1,u1,v1,w1,e1
    real(mytype), allocatable, dimension(:,:,:) :: r,u,v,w,e,pre,tem,mu,ka
    TYPE (DECOMP_INFO), intent(IN) :: part
    allocate(r(1-nHalo:part%xsz(1)+nHalo, 1-nHalo:part%xsz(2)+nHalo, 1-nHalo:part%xsz(3)+nHalo))
    allocate(u(1-nHalo:part%xsz(1)+nHalo, 1-nHalo:part%xsz(2)+nHalo, 1-nHalo:part%xsz(3)+nHalo))
    allocate(v(1-nHalo:part%xsz(1)+nHalo, 1-nHalo:part%xsz(2)+nHalo, 1-nHalo:part%xsz(3)+nHalo))
    allocate(w(1-nHalo:part%xsz(1)+nHalo, 1-nHalo:part%xsz(2)+nHalo, 1-nHalo:part%xsz(3)+nHalo))
    allocate(e(1-nHalo:part%xsz(1)+nHalo, 1-nHalo:part%xsz(2)+nHalo, 1-nHalo:part%xsz(3)+nHalo))
    allocate(pre(1-nHalo:part%xsz(1)+nHalo, 1-nHalo:part%xsz(2)+nHalo, 1-nHalo:part%xsz(3)+nHalo))
    allocate(tem(1-nHalo:part%xsz(1)+nHalo, 1-nHalo:part%xsz(2)+nHalo, 1-nHalo:part%xsz(3)+nHalo))
    allocate(mu(1-nHalo:part%xsz(1)+nHalo, 1-nHalo:part%xsz(2)+nHalo, 1-nHalo:part%xsz(3)+nHalo))
    allocate(ka(1-nHalo:part%xsz(1)+nHalo, 1-nHalo:part%xsz(2)+nHalo, 1-nHalo:part%xsz(3)+nHalo))
#if defined(BL)
  ! Blasius
  call initField_BL(part,x,z,r,u,v,w,e,pre,tem,mu,ka) 
#elif defined(CHA) 
  ! channel
  call initField_CHA(part,rho,u,v,w,e,pre,tem,mu,ka)
#elif defined(TGV)
  ! TGV
  call initField_TGV(part,rho,u,v,w,e,pre,tem,mu,ka)
#elif defined(1D) 
  ! 1D wave
  call initField_1D(part,rho,u,v,w,e,pre,tem,mu,ka)
#endif
    r1 = r(1:part%xsz(1), 1:part%xsz(2), 1:part%xsz(3))
    u1 = u(1:part%xsz(1), 1:part%xsz(2), 1:part%xsz(3))
    v1 = v(1:part%xsz(1), 1:part%xsz(2), 1:part%xsz(3))
    w1 = w(1:part%xsz(1), 1:part%xsz(2), 1:part%xsz(3))
    e1 = e(1:part%xsz(1), 1:part%xsz(2), 1:part%xsz(3))
    deallocate(r,u,v,w,e,pre,tem,mu,ka)
  end subroutine


! interpolate in all three directions
  subroutine interp3D(input, output, imax, jmax, kmax, inew, jnew, knew)
    use decomp_2d
    use mod_math
    implicit none
    integer :: imax,jmax,kmax,inew,jnew,knew,i,j,k, nHalo
    real(mytype) :: zinterp
    real(mytype), dimension(:,:,:) :: input, output
    real(mytype), allocatable, dimension(:,:,:) :: tmpArray1,tmpArray2, out_y,out_z
    real(mytype), allocatable, dimension(:) :: tmpInt1,tmpInt2
    TYPE (DECOMP_INFO) :: tmpPart
    
    allocate(out_y(part2%ysz(1),part2%ysz(2),part2%ysz(3)))
    allocate(out_z(part2%zsz(1),part2%zsz(2),part2%zsz(3)))
    call transpose_x_to_y(output, out_y, part2)
    call transpose_y_to_z(out_y,  out_z, part2)
    ! i-direction
    call decomp_info_init(inew,jmax,kmax,tmpPart)
    allocate(   tmpArray1(inew,tmpPart%xsz(2),tmpPart%xsz(3)))
    allocate(     tmpInt1(imax))
    allocate(     tmpInt2(imax))
    do k=1,tmpPart%xsz(3)
      do j=1,tmpPart%xsz(2)
        do i=1,imax
          tmpInt1(i) = input(i,j,k)
        enddo
        call spline(xold,tmpInt1,imax,tmpInt2)
        do i=1,inew
          call splint(xold,tmpInt1,tmpInt2,imax,xnew(i),tmpArray1(i,j,k))
        enddo
      enddo
    enddo
    allocate(tmpArray2(tmpPart%ysz(1),tmpPart%ysz(2),tmpPart%ysz(3)))
    call transpose_x_to_y(tmpArray1, tmpArray2, tmpPart)
    deallocate(tmpArray1)
    deallocate(  tmpInt1)
    deallocate(  tmpInt2)
    call decomp_info_finalize(tmpPart)
    ! j-direction
    call decomp_info_init(inew,jnew,kmax,tmpPart)
    allocate(tmpArray1(tmpPart%ysz(1),jnew,tmpPart%ysz(3)))
    allocate(tmpInt1(jmax))
    allocate(tmpInt2(jmax))
    if (jnew == 1) then 
      tmpArray1(:,1,:) = tmpArray2(:,1,:)
    else if (jmax /= 1) then
      do k=1,tmpPart%ysz(3)
        do i=1,tmpPart%ysz(1)
          do j=1,jmax
            tmpInt1(j) = tmpArray2(i,j,k)
          enddo
          call spline(yold,tmpInt1,jmax,tmpInt2)
          do j=1,jnew
           call splint(yold,tmpInt1,tmpInt2,jmax,ynew(j),tmpArray1(i,j,k))
          enddo
        enddo
      enddo
    else if (jmax == 1) then
      do k=1,tmpPart%ysz(3)
        do i=1,tmpPart%ysz(1)
          do j=1,jnew
            tmpArray1(i,j,k) = tmpArray2(i,jmax,k)
          enddo
        enddo
      enddo
    endif
    deallocate(tmpArray2)
    allocate(tmpArray2(tmpPart%zsz(1),tmpPart%zsz(2),tmpPart%zsz(3)))
    call transpose_y_to_z(tmpArray1, tmpArray2, tmpPart)
    deallocate(tmpArray1)
    deallocate(  tmpInt1)
    deallocate(  tmpInt2)
    call decomp_info_finalize(tmpPart)
    ! k-direction
    allocate(tmpInt1(kmax))
    allocate(tmpInt2(kmax))
    do j=1,part2%zsz(2)
      do i=1,part2%zsz(1)
        do k=1,kmax
          tmpInt1(k) = tmpArray2(i,j,k)
        enddo
        call spline(zold_global,tmpInt1,kmax,tmpInt2)
        do k=1,knew
          if (znew_global(k) < zold_global(kmax)) call splint(zold_global,tmpInt1,tmpInt2,kmax,znew_global(k),out_z(i,j,k))
        enddo
      enddo
    enddo
    deallocate(tmpArray2)
    deallocate(tmpInt1)
    deallocate(tmpInt2)
    call transpose_z_to_y(out_z, out_y, part2)
    call transpose_y_to_x(out_y, output,part2)
    deallocate(out_y)
    deallocate(out_z)
  end subroutine
end module
