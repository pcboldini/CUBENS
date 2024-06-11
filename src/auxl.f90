! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! auxiliary subroutines for parameters calculation and I/O
module mod_auxl
  use decomp_2d
  use mod_param
  use mod_grid
  use mod_finitediff
  implicit none
  contains
  
  ! calculation of the vorticity
  subroutine calcVort(vortx,vorty,vortz,strxz,u,v,w) 
    use decomp_2d
    implicit none
    integer i,j,k,c
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: vortx,vorty,vortz,strxz
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: u,v,w
    real(mytype) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, xa, za
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dwz,dwy,dwx,dvz,dvy,dvx,duz,duy,dux)
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          dux   = 0.0_mytype
          duy   = 0.0_mytype
          duz   = 0.0_mytype
          dvx   = 0.0_mytype
          dvy   = 0.0_mytype
          dvz   = 0.0_mytype
          dwx   = 0.0_mytype
          dwy   = 0.0_mytype
          dwz   = 0.0_mytype
          xa = xp(i)
          za = zp(k)
          !$acc loop seq
          do c = 1, nStencilVisc
            dux = dux + visc_ddx(c)*(u(i+c,j,k) - u(i-c,j,k))*xa
            duy = duy + visc_ddy(c)*(u(i,j+c,k) - u(i,j-c,k))
            duz = duz + visc_ddz(c)*(u(i,j,k+c) - u(i,j,k-c))*za
            dvx = dvx + visc_ddx(c)*(v(i+c,j,k) - v(i-c,j,k))*xa
            dvy = dvy + visc_ddy(c)*(v(i,j+c,k) - v(i,j-c,k))
            dvz = dvz + visc_ddz(c)*(v(i,j,k+c) - v(i,j,k-c))*za
            dwx = dwx + visc_ddx(c)*(w(i+c,j,k) - w(i-c,j,k))*xa
            dwy = dwy + visc_ddy(c)*(w(i,j+c,k) - w(i,j-c,k))
            dwz = dwz + visc_ddz(c)*(w(i,j,k+c) - w(i,j,k-c))*za
          enddo
          ! x: wall-normal, y: spanwise, z: streamwise
          vortx(i,j,k) = dwy - dvz ! du/dz - dw/dx -> vorticity wall-normal
          vorty(i,j,k) = duz - dwx ! dv/dx - du/dy -> vorticity spanwise
          vortz(i,j,k) = dvx - duy ! dw/dy - dv/dz -> vorticity streamwise
          strxz(i,j,k) = duz + dwx ! stress-tensor (2,1)=(1,2) 
        enddo
      enddo
    enddo
  end subroutine

! for channel flow: pressure gradient correction
  subroutine cmpbulkvel(r,w,wb)
    use decomp_2d
    use mod_param
    use mpi
    use mod_grid
    implicit none 
    integer ierr,i,j,k, nxnynz
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: w,r
    real(mytype) wb, da
    ! initialization
    wb  = 0.0_mytype
    !$acc parallel default(present) copy(wb)
    !$acc loop gang, vector collapse(3) reduction(+:wb) 
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=2,xsize(1)
             da = x(i)-x(i-1)
             wb = wb  + 0.5_mytype*(r(i,j,k)*w(i,j,k)+r(i-1,j,k)*w(i-1,j,k))*da
          enddo
       enddo
    enddo
    !$acc end parallel
    nxnynz = jmax*kmax*len_x
    call mpi_allreduce(MPI_IN_PLACE, wb, 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   wb = wb/nxnynz
  end subroutine cmpbulkvel

! calculation of bulk properties for quick monitoring of the simulation
  subroutine cmpbulk(istep,wt1,time,dt,CFL_new,rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz,wb,dp)
    use decomp_2d
    use mod_param
    use mpi
    use mod_grid
    use mod_eos
    implicit none 
    integer ierr,i,j,k,istep,ioutput,nxnynz
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz
    real(mytype) ddx,ddy,ddz,time,dt,isNan,isNanGlobal,CFL_new
    real(mytype) rb,wb,preb,dp,kib,eb,str,tmp,enst,sos,max_Mach
    real(8) :: wt1
    ! initialization
    isNan = 0.0_mytype
    isNanGlobal = 0.0_mytype
    rb   = 0.0_mytype
    wb   = 0.0_mytype
    preb = 0.0_mytype
    kib  = 0.0_mytype
    eb   = 0.0_mytype
    str  = 0.0_mytype
    enst = 0.0_mytype
    max_Mach = 0.0_mytype
    ! domain volume
    nxnynz = len_x*len_z*len_y
    ! grid spacing in y- and z-direction
    if (xsize(2) == 1) then
      ddy = 1.0_mytype
    else
      ddy = y(2)-y(1)
    endif
    ddz = z(2)-z(1)
    ! calculation of bulk properties
    !$acc parallel default(present) copy(rb,wb,preb,kib,eb,enst,max_Mach)
    !$acc loop gang, vector collapse(3) reduction(+:rb,wb,preb,kib,eb,enst) reduction(max:max_Mach)
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=2,xsize(1)
          ddx = x(i)-x(i-1)
          ! density
          rb  = rb  + 0.5_mytype*(rho(i,j,k)+rho(i-1,j,k))*ddx
          ! streamwise velocity
          wb  = wb  + 0.5_mytype*(w(i,j,k)+w(i-1,j,k))*ddx
          ! pressure
          preb = preb + 0.5_mytype*(pre(i,j,k)+pre(i-1,j,k))*ddx
          ! kinetic energy
          kib = kib + 0.5_mytype*rho(i,j,k)*(w(i,j,k)**2+v(i,j,k)**2+u(i,j,k)**2)*ddx
          ! internal energy
          eb  = eb  + 0.5_mytype*(ien(i,j,k)+ien(i-1,j,k))**2*ddx
          ! enstrophy
          enst = enst + mu(i,j,k)*(vortx(i,j,k)**2 + vorty(i,j,k)**2 + vortz(i,j,k)**2)*ddx
          call calcSOS(rho(i,j,k),ien(i,j,k),sos)
          ! Mach number
          max_Mach = max(sqrt(w(i,j,k)**2+v(i,j,k)**2+u(i,j,k)**2)/sos, max_Mach)
        enddo
      enddo
    enddo
    !$acc end parallel
    ! Wall shear stress
    !$acc parallel default(present) copy(str)
    !$acc loop gang, vector reduction(+:str)
    do k=1,xsize(3)
      do j=1,xsize(2)
        ! approximation
        str = str + (mu(2,j,k)*w(2,j,k)/x(2)) 
      enddo
    enddo
    !$acc end parallel
    rb = rb*ddy*ddz
    wb = wb*ddy*ddz
    preb = preb*ddy*ddz
    kib = kib*ddy*ddz
    eb = eb*ddy*ddz
    enst = enst*ddy*ddz
    ! combine all MPI processes
    call mpi_allreduce(MPI_IN_PLACE, rb  , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   rb   = rb/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, wb  , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   wb   = wb/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, preb  , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr); preb = preb/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, kib , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   kib  = kib/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, eb  , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   eb   = eb/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, str , 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   str  = str/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, enst, 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr);   enst = enst/nxnynz
    call mpi_allreduce(MPI_IN_PLACE, max_Mach, 1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr);   
    ! calculate maximal Mach number
    ! print the step, timestep, CFL, and the previously calculated bulk quantities
    if  (nrank .eq. 0) then
#if defined(CHA)      
      if ((mod(istep,20*intvPrint).eq.0)) then 
        write(*,'(A7,14A15)') 'step', 'wall-time', 'dt',  'CFL', 'time', 'wall-stress', &
                              'bulk w-vel', 'bulk rho', 'bulk P', 'dpdz', 'ke-bulk', 'eint-bulk' , 'enstrophy', 'max_Mach'
      endif
      if (mod(istep,intvPrint).eq.0) write(*,'(I7,13E15.6)') &
                     istep, MPI_WTIME()-wt1, dt, CFL_new, time, str, wb, rb, preb, dp, kib, eb, enst, max_Mach
#else
      if ((mod(istep,20*intvPrint).eq.0)) then 
        write(*,'(A7,13A16)') 'step', 'wall-time', 'dt',  'CFL', 'time', 'wall-stress', &
                              'bulk w-vel', 'bulk rho', 'bulk P', 'ke-bulk', 'eint-bulk' , 'enstrophy', 'max_Mach'
      endif
      if (mod(istep,intvPrint).eq.0) write(*,'(I7,13E16.7)') &
                     istep, MPI_WTIME()-wt1, dt, CFL_new, time, str, wb, rb, preb, kib, eb, enst, max_Mach
#endif
    endif
    ! check if there is any NaN, in case stop the simulation
    if (wb .ne. wb) isNan = 1
    call mpi_allreduce(isNan, isNanGlobal,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (isNanGlobal>0) then
      call decomp_2d_finalize
      call mpi_finalize(ierr)
      stop
    endif
  end subroutine

  ! calculate statistics
  subroutine calcStats(qave,factAvg,countAvg,rho,u,v,w,ien,pre,tem,mu,ka)
    use decomp_2d
    use mod_param
    implicit none
    integer i,j,k,h,c,countAvg
    real(mytype) :: factAvg
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
    real(mytype) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, xa, za
    real(mytype) :: dilla2, dtx, dty, dtz
#if defined(BL) || defined(TGV)
    real(mytype), dimension(:,:,:), intent(inout) :: qave
    real(mytype) :: rho_2d,u_2d,v_2d,w_2d,pre_2d,tem_2d,ien_2d,mu_2d,ka_2d
    real(mytype) :: tauxx_2d,tauxy_2d,tauxz_2d,tauyy_2d,tauyz_2d,tauzz_2d,qx_2d,qy_2d,qz_2d
#elif defined(CHA) 
    real(mytype), dimension(:,:), intent(inout) :: qave
    real(mytype) :: rho_1d,u_1d,v_1d,w_1d,pre_1d,tem_1d,ien_1d,mu_1d,ka_1d
    real(mytype) :: tauxx_1d,tauxy_1d,tauxz_1d,tauyy_1d,tauyz_1d,tauzz_1d,qx_1d,qy_1d,qz_1d
#endif
#if defined(BL)
    ! for boundary layer
    !$acc parallel loop gang collapse(2) default(present)
    do k=1,xsize(3)
      do i=1,xsize(1)
        rho_2d = 0.0_mytype
        u_2d = 0.0_mytype
        v_2d = 0.0_mytype
        w_2d = 0.0_mytype
        pre_2d = 0.0_mytype
        tem_2d = 0.0_mytype
        ien_2d = 0.0_mytype
        mu_2d = 0.0_mytype
        ka_2d = 0.0_mytype
        tauxx_2d = 0.0_mytype
        tauxy_2d = 0.0_mytype
        tauxz_2d = 0.0_mytype
        tauyy_2d = 0.0_mytype
        tauyz_2d = 0.0_mytype
        tauzz_2d = 0.0_mytype
        qx_2d = 0.0_mytype
        qy_2d = 0.0_mytype
        qz_2d = 0.0_mytype
        !$acc loop reduction(+:rho_2d,u_2d,v_2d,w_2d,pre_2d,tem_2d,ien_2d,mu_2d,ka_2d) &
        !$acc reduction(+:tauxx_2d,tauxy_2d,tauxz_2d,tauyy_2d,tauyz_2d,tauzz_2d,qx_2d,qy_2d,qz_2d) &
        !$acc private(dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz,dtx,dty,dtz) 
        do j=1,xsize(2)
          rho_2d = rho_2d + rho(i,j,k)
          u_2d = u_2d + u(i,j,k)
          v_2d = v_2d + v(i,j,k)
          w_2d = w_2d + w(i,j,k)
          pre_2d = pre_2d + pre(i,j,k)
          tem_2d = tem_2d + tem(i,j,k)
          ien_2d = ien_2d + ien(i,j,k)
          mu_2d = mu_2d + mu(i,j,k)
          ka_2d = ka_2d + ka(i,j,k)
          dux   = 0.0_mytype
          duy   = 0.0_mytype
          duz   = 0.0_mytype
          dvx   = 0.0_mytype
          dvy   = 0.0_mytype
          dvz   = 0.0_mytype
          dwx   = 0.0_mytype
          dwy   = 0.0_mytype
          dwz   = 0.0_mytype
          dtx   = 0.0_mytype
          dty   = 0.0_mytype
          dtz   = 0.0_mytype
          xa = xp(i)
          za = zp(k)
          !$acc loop seq
          do c = 1, nStencilVisc
            dux = dux + visc_ddx(c)*(u(i+c,j,k) - u(i-c,j,k))*xa
            duy = duy + visc_ddy(c)*(u(i,j+c,k) - u(i,j-c,k))
            duz = duz + visc_ddz(c)*(u(i,j,k+c) - u(i,j,k-c))*za
            dvx = dvx + visc_ddx(c)*(v(i+c,j,k) - v(i-c,j,k))*xa
            dvy = dvy + visc_ddy(c)*(v(i,j+c,k) - v(i,j-c,k))
            dvz = dvz + visc_ddz(c)*(v(i,j,k+c) - v(i,j,k-c))*za
            dwx = dwx + visc_ddx(c)*(w(i+c,j,k) - w(i-c,j,k))*xa
            dwy = dwy + visc_ddy(c)*(w(i,j+c,k) - w(i,j-c,k))
            dwz = dwz + visc_ddz(c)*(w(i,j,k+c) - w(i,j,k-c))*za
            dtx = dtx + visc_ddx(c)*(tem(i+c,j,k) - tem(i-c,j,k))*xa
            dty = dty + visc_ddy(c)*(tem(i,j+c,k) - tem(i,j-c,k))
            dtz = dtz + visc_ddz(c)*(tem(i,j,k+c) - tem(i,j,k-c))*za
          enddo
          dilla2   = dux + dvy + dwz
          tauxx_2d = tauxx_2d + mu(i,j,k)*(2.0_mytype*dux - 2.0_mytype/3.0_mytype*dilla2) 
          tauxy_2d = tauxy_2d + mu(i,j,k)*(duy + dvx)
          tauxz_2d = tauxz_2d + mu(i,j,k)*(duz + dwx)
          tauyy_2d = tauyy_2d + mu(i,j,k)*(2.0_mytype*dvy - 2.0_mytype/3.0_mytype*dilla2)
          tauyz_2d = tauyz_2d + mu(i,j,k)*(dvz + dwy)
          tauzz_2d = tauzz_2d + mu(i,j,k)*(2.0_mytype*dwz - 2.0_mytype/3.0_mytype*dilla2)
          qx_2d    = qx_2d    + ka(i,j,k)*dtx
          qy_2d    = qy_2d    + ka(i,j,k)*dty
          qz_2d    = qz_2d    + ka(i,j,k)*dtz
        enddo
      ! primary variables
      qave(i,k,1) = rho_2d  
      qave(i,k,2) = u_2d  
      qave(i,k,3) = v_2d 
      qave(i,k,4) = w_2d
      qave(i,k,5) = pre_2d
      qave(i,k,6) = tem_2d
      qave(i,k,7) = ien_2d
      qave(i,k,8) = mu_2d
      qave(i,k,9) = ka_2d
      ! double products
      qave(i,k,10) = rho_2d*u_2d
      qave(i,k,11) = rho_2d*v_2d
      qave(i,k,12) = rho_2d*w_2d
      qave(i,k,13) = rho_2d*tem_2d
      qave(i,k,14) = u_2d*u_2d
      qave(i,k,15) = u_2d*v_2d
      qave(i,k,16) = u_2d*w_2d
      qave(i,k,17) = v_2d*v_2d
      qave(i,k,18) = v_2d*w_2d
      qave(i,k,19) = w_2d*w_2d
      ! triple products
      qave(i,k,20) = rho_2d*u_2d*u_2d
      qave(i,k,21) = rho_2d*u_2d*v_2d
      qave(i,k,22) = rho_2d*u_2d*w_2d
      qave(i,k,23) = rho_2d*v_2d*v_2d
      qave(i,k,24) = rho_2d*v_2d*w_2d
      qave(i,k,25) = rho_2d*w_2d*w_2d
      ! stress tensor
      qave(i,k,26) = tauxx_2d
      qave(i,k,27) = tauxy_2d 
      qave(i,k,28) = tauxz_2d 
      qave(i,k,29) = tauyy_2d 
      qave(i,k,30) = tauyz_2d 
      qave(i,k,31) = tauzz_2d
      ! heat flux 
      qave(i,k,32) = qx_2d
      qave(i,k,33) = qy_2d
      qave(i,k,34) = qz_2d
    enddo
    enddo
    factAvg = factAvg + 1.0_mytype*xsize(2)
#elif defined(CHA)
    ! for channel
    !$acc parallel loop gang default(present)
    do i=1,xsize(1)
      rho_1d = 0.0_mytype
      u_1d = 0.0_mytype
      v_1d = 0.0_mytype
      w_1d = 0.0_mytype
      pre_1d = 0.0_mytype
      tem_1d = 0.0_mytype
      ien_1d = 0.0_mytype
      mu_1d = 0.0_mytype
      ka_1d = 0.0_mytype
      tauxx_1d = 0.0_mytype
      tauxy_1d = 0.0_mytype
      tauxz_1d = 0.0_mytype
      tauyy_1d = 0.0_mytype
      tauyz_1d = 0.0_mytype
      tauzz_1d = 0.0_mytype
      qx_1d = 0.0_mytype
      qy_1d = 0.0_mytype
      qz_1d = 0.0_mytype
      !$acc loop collapse(2) reduction(+:rho_1d,u_1d,v_1d,w_1d,pre_1d,tem_1d,ien_1d,mu_1d,ka_1d) &
      !$acc reduction(+:tauxx_1d,tauxy_1d,tauxz_1d,tauyy_1d,tauyz_1d,tauzz_1d,qx_1d,qy_1d,qz_1d) &
      !$acc private(dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz,dtx,dty,dtz) 
      do k=1,xsize(3)
        do j=1,xsize(2)
          rho_1d = rho_1d + rho(i,j,k)
          u_1d = u_1d + u(i,j,k)
          v_1d = v_1d + v(i,j,k)
          w_1d = w_1d + w(i,j,k)
          pre_1d = pre_1d + pre(i,j,k)
          tem_1d = tem_1d + tem(i,j,k)
          ien_1d = ien_1d + ien(i,j,k)
          mu_1d = mu_1d + mu(i,j,k)
          ka_1d = ka_1d + ka(i,j,k)
          dux   = 0.0_mytype
          duy   = 0.0_mytype
          duz   = 0.0_mytype
          dvx   = 0.0_mytype
          dvy   = 0.0_mytype
          dvz   = 0.0_mytype
          dwx   = 0.0_mytype
          dwy   = 0.0_mytype
          dwz   = 0.0_mytype
          dtx   = 0.0_mytype
          dty   = 0.0_mytype
          dtz   = 0.0_mytype
          xa = xp(i)
          za = zp(k)
          !$acc loop seq
          do c = 1, nStencilVisc
            dux = dux + visc_ddx(c)*(u(i+c,j,k) - u(i-c,j,k))*xa
            duy = duy + visc_ddy(c)*(u(i,j+c,k) - u(i,j-c,k))
            duz = duz + visc_ddz(c)*(u(i,j,k+c) - u(i,j,k-c))*za
            dvx = dvx + visc_ddx(c)*(v(i+c,j,k) - v(i-c,j,k))*xa
            dvy = dvy + visc_ddy(c)*(v(i,j+c,k) - v(i,j-c,k))
            dvz = dvz + visc_ddz(c)*(v(i,j,k+c) - v(i,j,k-c))*za
            dwx = dwx + visc_ddx(c)*(w(i+c,j,k) - w(i-c,j,k))*xa
            dwy = dwy + visc_ddy(c)*(w(i,j+c,k) - w(i,j-c,k))
            dwz = dwz + visc_ddz(c)*(w(i,j,k+c) - w(i,j,k-c))*za
            dtx = dtx + visc_ddx(c)*(tem(i+c,j,k) - tem(i-c,j,k))*xa
            dty = dty + visc_ddy(c)*(tem(i,j+c,k) - tem(i,j-c,k))
            dtz = dtz + visc_ddz(c)*(tem(i,j,k+c) - tem(i,j,k-c))*za
          enddo
          dilla2   = dux + dvy + dwz
          tauxx_1d = tauxx_1d + mu(i,j,k)*(2.0_mytype*dux - 2.0_mytype/3.0_mytype*dilla2) 
          tauxy_1d = tauxy_1d + mu(i,j,k)*(duy + dvx)
          tauxz_1d = tauxz_1d + mu(i,j,k)*(duz + dwx)
          tauyy_1d = tauyy_1d + mu(i,j,k)*(2.0_mytype*dvy - 2.0_mytype/3.0_mytype*dilla2)
          tauyz_1d = tauyz_1d + mu(i,j,k)*(dvz + dwy)
          tauzz_1d = tauzz_1d + mu(i,j,k)*(2.0_mytype*dwz - 2.0_mytype/3.0_mytype*dilla2)
          qx_1d    = qx_1d    + ka(i,j,k)*dtx
          qy_1d    = qy_1d    + ka(i,j,k)*dty
          qz_1d    = qz_1d    + ka(i,j,k)*dtz
        enddo
      enddo
      ! primary variables
      qave(i,1) = rho_1d  
      qave(i,2) = u_1d  
      qave(i,3) = v_1d 
      qave(i,4) = w_1d
      qave(i,5) = pre_1d
      qave(i,6) = tem_1d
      qave(i,7) = ien_1d
      qave(i,8) = mu_1d
      qave(i,9) = ka_1d
      ! double products
      qave(i,10) = rho_1d*u_1d
      qave(i,11) = rho_1d*v_1d
      qave(i,12) = rho_1d*w_1d
      qave(i,13) = rho_1d*tem_1d
      qave(i,14) = u_1d*u_1d
      qave(i,15) = u_1d*v_1d
      qave(i,16) = u_1d*w_1d
      qave(i,17) = v_1d*v_1d
      qave(i,18) = v_1d*w_1d
      qave(i,19) = w_1d*w_1d
      ! triple products
      qave(i,20) = rho_1d*u_1d*u_1d
      qave(i,21) = rho_1d*u_1d*v_1d
      qave(i,22) = rho_1d*u_1d*w_1d
      qave(i,23) = rho_1d*v_1d*v_1d
      qave(i,24) = rho_1d*v_1d*w_1d
      qave(i,25) = rho_1d*w_1d*w_1d
      ! stress tensor
      qave(i,26) = tauxx_1d
      qave(i,27) = tauxy_1d 
      qave(i,28) = tauxz_1d 
      qave(i,29) = tauyy_1d 
      qave(i,30) = tauyz_1d 
      qave(i,31) = tauzz_1d
      ! heat flux 
      qave(i,32) = qx_1d
      qave(i,33) = qy_1d
      qave(i,34) = qz_1d
    enddo
    factAvg = factAvg + 1.0_mytype*xsize(2)*xsize(3)
#endif  
    countAvg = countAvg + 1
  end subroutine

! I/O planes
  subroutine output2dPlane(part,nHaloIn,istep,dir,loc,tmp01,name011,name012,tmp02,name021,name022,tmp03,name031,name032,tmp04,   &
                                                      name041,name042,tmp05,name051,name052,tmp06,name061,name062,tmp07,name071, &
                                                      name072,tmp08,name081,name082,tmp09,name091,name092,tmp10,name101,name102, &
                                                      tmp11,name111,name112,tmp12,name121,name122)
    ! -- dir specifies the direction of the plane:
    ! - 1 --> x-normal plane
    ! - 2 --> y-normal plane
    ! - 3 --> z-normal plane
    ! loc is the location of the plane in global coordinates
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,dir,loc,nHaloIn
    character*7 cha
    character(len=1024) :: cha2
    character(len=*), optional :: name011,name021,name031,name041,name051,name061,name071,name081,name091,name101,name111,name121
    character(len=*), optional :: name012,name022,name032,name042,name052,name062,name072,name082,name092,name102,name112,name122
    real(mytype), dimension(1-nHaloIn:,1-nHaloIn:,1-nHaloIn:), optional :: &
                        tmp01,tmp02,tmp03,tmp04,tmp05,tmp06,tmp07,tmp08,tmp09,tmp10,tmp11,tmp12
    TYPE (DECOMP_INFO) :: part                    
    real(mytype), allocatable, dimension(:,:,:) :: tmp
    allocate(tmp(part%xsz(1),part%xsz(2),part%xsz(3)))
    write(cha,'(I0.7)') istep
    write(cha2,'(I0)') loc
    ! all the planes are written in ./output/planes/
    if(present(tmp01)) then 
      tmp = tmp01(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name011)//trim(cha2)//'.'//trim(name012)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp02)) then 
      tmp = tmp02(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))  
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name021)//trim(cha2)//'.'//trim(name022)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp03)) then 
      tmp = tmp03(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name031)//trim(cha2)//'.'//trim(name032)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp04)) then 
      tmp = tmp04(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) 
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name041)//trim(cha2)//'.'//trim(name042)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp05)) then 
      tmp = tmp05(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name051)//trim(cha2)//'.'//trim(name052)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp06)) then 
      tmp =  tmp06(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) 
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name061)//trim(cha2)//'.'//trim(name062)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp07)) then 
      tmp =  tmp07(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) 
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name071)//trim(cha2)//'.'//trim(name072)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp08)) then 
      tmp =  tmp08(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name081)//trim(cha2)//'.'//trim(name082)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp09)) then 
      tmp =  tmp09(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name091)//trim(cha2)//'.'//trim(name092)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp10)) then 
      tmp =  tmp10(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name101)//trim(cha2)//'.'//trim(name102)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp11)) then 
      tmp =  tmp11(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name111)//trim(cha2)//'.'//trim(name112)//cha//'.bin', &
                                 'dummy',part)
    endif
    if(present(tmp12)) then 
      tmp =  tmp12(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3))
      call decomp_2d_write_plane(1,tmp,dir,loc,'.','output/planes/'//trim(name121)//trim(cha2)//'.'//trim(name122)//cha//'.bin', &
                                 'dummy',part)
    endif
    if (nrank == 0) then
      if (istep == 0) then
        if (dir==1) then 
          write(stdout,* ) 'Initial x-plane:'
        else if (dir==2) then 
          write(stdout,* ) 'Initial y-plane:'
        else if (dir==3) then
          write(stdout,* ) 'Initial z-plane:'
        endif
        write(stdout,'(A)') 'o--------------------------------------------------o'
        write(stdout,'(A)') 'Calculation and writing                    done!'
        write(stdout,* ) 
      else if ( (istep > 0) .AND. (dir==1)) then 
        write(stdout,'(A, I10, A, I5)') 'Writing x-plane at step: ', istep,' and index:', loc
      else if ( (istep > 0) .AND. (dir==2)) then 
        write(stdout,'(A, I10, A, I5)') 'Writing y-plane at step: ', istep,' and index:', loc
      else if ( (istep > 0) .AND. (dir==3)) then
        write(stdout,'(A, I10, A, I5)') 'Writing z-plane at step: ', istep,' and index:', loc
      endif
    endif
    deallocate(tmp)
  end subroutine

! I/O load restart
  subroutine loadRestart(istep,time,rho,u,v,w,ien,nHaloIn,part)
    use mpi
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,ierr,nHaloIn
    character*7 cha
    real(mytype), dimension(1-nHaloIn:,1-nHaloIn:,1-nHaloIn:) :: rho,u,v,w,ien
    real(mytype), allocatable, dimension(:,:,:) :: tmp
    real(mytype), dimension(5) :: infoRestart
    real(mytype) :: time
    LOGICAL :: file_exists
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    TYPE(DECOMP_INFO), intent(IN) :: part
    file_exists = .FALSE.
    allocate(tmp(part%xsz(1),part%xsz(2),part%xsz(3)))
    write(cha,'(I0.7)') istep
    ! all the restart are located in ./output/restart
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'output/restart/ruvwe.'//cha//'.bin', &
                       MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
    call MPI_FILE_GET_SIZE(fh,filesize,ierr)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_read_scalar(fh,disp,5,infoRestart)
    INQUIRE (FILE='output/restart/ruvwe.'//cha//'.bin', EXIST=file_exists)
    ! check if restart file exists
    if (file_exists) then
      if(nrank==0) write(stdout,'(A)') 'restart FILE exists!'
    else
      if(nrank==0) write(stdout,'(A)') '!restart FILE does not exist!'
      call MPI_FILE_CLOSE(fh,ierr)
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      stop
    end if
    ! check if the restart file precision is correct
    if (infoRestart(1).ne.1.0*mytype) then 
      if (nrank.eq.0) write(stdout,'(A)') 'Simulation stopped!'
      if (nrank.eq.0) write(stdout,'(A)') 'Restart file has different precision.'
      call MPI_FILE_CLOSE(fh,ierr)
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      stop
    endif
    ! check if the restart file has the correct size
    if ((infoRestart(2).ne.1.0*nx_global) .or. &
        (infoRestart(3).ne.1.0*ny_global) .or. &
        (infoRestart(4).ne.1.0*nz_global)) then 
      if (nrank.eq.0) write(stdout,*) 'Simulation stopped!'
      if (nrank.eq.0) write(stdout,*) 'Restart has imax/jmax/kmax : ', nint(infoRestart(2)), &
                                                                       nint(infoRestart(3)), &
                                                                       nint(infoRestart(4))
      if (nrank.eq.0) write(stdout,*) 'This sim has imax/jmax/kmax: ', nx_global,ny_global,nz_global
      call MPI_FILE_CLOSE(fh,ierr)
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      stop
    endif 
    ! open the five primary non-conservative variables: density, velocities, and internal energy
    call decomp_2d_read_var(fh,disp,1,tmp,part); rho(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call decomp_2d_read_var(fh,disp,1,tmp,part);   u(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call decomp_2d_read_var(fh,disp,1,tmp,part);   v(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call decomp_2d_read_var(fh,disp,1,tmp,part);   w(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call decomp_2d_read_var(fh,disp,1,tmp,part); ien(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3)) = tmp
    call MPI_FILE_CLOSE(fh,ierr)
    deallocate(tmp)
    time = infoRestart(5)
    if(nrank==0) write(stdout,'(A,E12.5,A,I0)') 'The simulation restarts at t: ',time, ' and at step: ', istep
    if(nrank==0) write(stdout,*) ''
  end subroutine

! I/O save restart
  subroutine saveRestart(istep,time,rho,u,v,w,ien,nHaloIn,part,interpol)
    use mpi
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,ierr,nHaloIn
    character*7 cha
    real(mytype), dimension(1-nHaloIn:,1-nHaloIn:,1-nHaloIn:) :: rho,u,v,w,ien
    real(mytype), allocatable, dimension(:,:,:) :: tmp
    real(mytype), dimension(5) :: infoRestart
    real(mytype) :: time
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    TYPE(DECOMP_INFO), intent(IN) :: part
    character(len=*), optional :: interpol
    ! read the five infoRestart arguments
    infoRestart(1) = 1.0*mytype
    ! distinguish interpolated file
    if(interpol == 'aI') then 
      infoRestart(2) = 1.0*inew
      infoRestart(3) = 1.0*jnew
      infoRestart(4) = 1.0*knew
    else
      infoRestart(2) = 1.0*imax
      infoRestart(3) = 1.0*jmax
      infoRestart(4) = 1.0*kmax
    endif
    infoRestart(5) = time 
    allocate(tmp(part%xsz(1),part%xsz(2),part%xsz(3)))
    write(cha,'(I0.7)') istep
    if (nrank == 0) then
      if (istep == 0) then
        write(stdout,* ) 
        write(stdout,* ) 'Initial restart condition'
        write(stdout,'(A)') 'o--------------------------------------------------o'
        write(stdout,'(A)') 'Calculation and writing                    done!'
        write(stdout,* ) 
      else if ( istep > 0 ) then 
        write(stdout,'(A, I10)') 'Writing restart at step:              ', istep    
      endif
    endif
    ! all the restart files are written in ./restart
    if(interpol == 'bI') then 
      if (nrank == 0) write(stdout,* ) 'before interpolation'
      call MPI_FILE_OPEN(MPI_COMM_WORLD, 'output/restart/ruvwe.'//trim(interpol)//'.'//cha//'.bin', &
                        MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    else
      call MPI_FILE_OPEN(MPI_COMM_WORLD, 'output/restart/ruvwe.'//cha//'.bin', &
                        MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    endif
    ! writing of the restart file
    filesize = 0_MPI_OFFSET_KIND
    ! guarantee overwriting
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_scalar(fh,disp,5,infoRestart)
    ! write the five primary non-conservative variables: density, velocities, and internal energy
    tmp = rho(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    tmp =   u(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    tmp =   v(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    tmp =   w(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    tmp = ien(1:part%xsz(1),1:part%xsz(2),1:part%xsz(3));  call decomp_2d_write_var(fh,disp,1,tmp,part)
    call MPI_FILE_CLOSE(fh,ierr)
    deallocate(tmp)
  end subroutine
  ! I/O save statistics
  subroutine saveStats(part,istep,dt,qave,factAvg,countAvg)
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,countAvg
    character*7 cha
    real(mytype) :: dt, factAvg, factAvg_inv
    TYPE (DECOMP_INFO) :: part
    logical :: exist
#if defined(BL) || defined(TGV)
    ! for boundary layer
    real(mytype), allocatable, dimension(:,:,:) :: tmpPlane
    real(mytype), dimension(:,:,:), intent(in) :: qave
    allocate(tmpPlane(part%xsz(1),part%xsz(2),part%xsz(3)))
    tmpPlane = 0.0_mytype
#elif defined(CHA)
    real(mytype), allocatable, dimension(:,:) :: tmpCut 
    real(mytype), dimension(:,:), intent(in) :: qave
#endif
    write(cha,'(I0.7)') istep
    factAvg_inv = 1.0_mytype/factAvg
    ! write stats planes in /output/stats/ 
    ! for boundary layer
#if defined(BL)
    ! density (1)
    tmpPlane(:,1,:) = qave(:,:,1)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/r_avg.'//cha//'.bin','dummy')
    ! u-velocity (2)
    tmpPlane(:,1,:) = qave(:,:,2)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/u_avg.'//cha//'.bin','dummy')
    ! v-velocity (3)
    tmpPlane(:,1,:) = qave(:,:,3)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/v_avg.'//cha//'.bin','dummy')
    ! w-velocity (4) 
    tmpPlane(:,1,:) = qave(:,:,4)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/w_avg.'//cha//'.bin','dummy')
    ! pressure (5)
    tmpPlane(:,1,:) = qave(:,:,5)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/p_avg.'//cha//'.bin','dummy')
    ! temperature (6)   
    tmpPlane(:,1,:) = qave(:,:,6)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/t_avg.'//cha//'.bin','dummy')
    ! viscosity (8)
    tmpPlane(:,1,:) = qave(:,:,8)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/mu_avg.'//cha//'.bin','dummy') 
    ! conductivity (9)
    tmpPlane(:,1,:) = qave(:,:,9)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/ka_avg.'//cha//'.bin','dummy') 
    ! ru (10)
    tmpPlane(:,1,:) = qave(:,:,10)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/ru_avg.'//cha//'.bin','dummy') 
    ! rv (11) 
    tmpPlane(:,1,:) = qave(:,:,11)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/rv_avg.'//cha//'.bin','dummy') 
    ! rw (12)
    tmpPlane(:,1,:) = qave(:,:,12)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/rw_avg.'//cha//'.bin','dummy') 
    ! rt (13)
    tmpPlane(:,1,:) = qave(:,:,13)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/rt_avg.'//cha//'.bin','dummy')   
    ! ruu (20)
    tmpPlane(:,1,:) = qave(:,:,20)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/ruu_avg.'//cha//'.bin','dummy')   
    ! ruv (21)
    tmpPlane(:,1,:) = qave(:,:,21)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/ruv_avg.'//cha//'.bin','dummy')   
    ! ruw (22)
    tmpPlane(:,1,:) = qave(:,:,22)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/ruw_avg.'//cha//'.bin','dummy')   
    ! rvv (23)
    tmpPlane(:,1,:) = qave(:,:,23)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/rvv_avg.'//cha//'.bin','dummy')    
    ! rvw (24)
    tmpPlane(:,1,:) = qave(:,:,24)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/rvw_avg.'//cha//'.bin','dummy')    
    ! rww (25)
    tmpPlane(:,1,:) = qave(:,:,25)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/rww_avg.'//cha//'.bin','dummy')    
    ! tau_xx (26)
    tmpPlane(:,1,:) = qave(:,:,26)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/tauxx_avg.'//cha//'.bin','dummy')    
    ! tau_xy (27)
    tmpPlane(:,1,:) = qave(:,:,27)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/tauxy_avg.'//cha//'.bin','dummy')    
    ! tau_xz (28) 
    tmpPlane(:,1,:) = qave(:,:,28)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/tauxz_avg.'//cha//'.bin','dummy')    
    ! tau_yy (29)
    tmpPlane(:,1,:) = qave(:,:,29)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/tauyy_avg.'//cha//'.bin','dummy')    
    ! tau_yz (30)
    tmpPlane(:,1,:) = qave(:,:,30)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/tauyz_avg.'//cha//'.bin','dummy')    
    ! tau_zz (31) 
    tmpPlane(:,1,:) = qave(:,:,31)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/tauzz_avg.'//cha//'.bin','dummy')
    ! q_x (32)
    tmpPlane(:,1,:) = qave(:,:,32)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/qx_avg.'//cha//'.bin','dummy')
    ! q_y (33)
    tmpPlane(:,1,:) = qave(:,:,33)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/qy_avg.'//cha//'.bin','dummy')
    ! q_z (34)
    tmpPlane(:,1,:) = qave(:,:,34)*factAvg_inv; &
    call decomp_2d_write_plane(1,tmpPlane,2,1, '.','output/stats/qz_avg.'//cha//'.bin','dummy')
    ! for channel
#elif defined(CHA) 
    ! call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(1),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !   if(myid == 0) then
    !     open(newunit=iunit,file=fname)
    !     do i=1,ng(1)
    !       write(iunit,fmt_rp) (i-.5)*dl(1),p1d(i)
    !     end do
    !     close(iunit)
    !   end if
    ! tmpPlane(:,1,:) = qave(:,:,1)*factAvg_inv; &
    ! call decomp_2d_write_plane(1,tmpPlane,2,1, '.','postproc/stats/r_avg.'//cha//'.bin','dummy')
#endif        
    ! save information for this average
    if (nrank == 0) then
#if defined(BL)
      write(stdout,'(A, I0)') 'Boundary layer: writing time- and span-averaged statistics at final step: ', istep
#elif defined(CHA)
      write(stdout,'(A, I0)') 'Channel: writing time-, span- and stream-averaged statistics at final step: ', istep
#endif          
      inquire(file="postproc/stats_info.txt", exist=exist)
      if (exist) then
        open(11,file = 'postproc/stats_info.txt',status="old",position="append",action="write")
      else
        open(11,file = 'postproc/stats_info.txt',status="new", action="write")
      endif   
      write(11,'(A, I0)') 'step: ', istep
      write(11,'(A, I0)') 'number of averaged files: ', countAvg
      write(11,'(A, I0)') 'average factor: ', int(factAvg)
      write(11,'(A, E17.6)') 'time from begin averaging: ', dt*(istep - saveStatsAfter)
      close(11)
    endif
  end subroutine
end module
