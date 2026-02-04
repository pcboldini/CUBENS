! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! auxiliary subroutines for parameters calculation and I/O

module mod_auxl
  use decomp_2d
  use mod_param
  use mod_grid
  use mod_finitediff
  use mod_eos
  use mod_eos_var
  use mod_eos_visc
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
    integer ierr,i,j,k
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: w,r
    real(mytype) wb, da, nxnynz
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
    character(len=100) :: filename
    integer ierr,i,j,k,istep,ioutput
    real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka,vortx,vorty,vortz
    real(mytype) ddx,ddy,ddz,time,dt,isNan,isNanGlobal,CFL_new,nxnynz
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
#if defined(TGV)
    ! Specify the filename
    filename = 'tests/TGV/output.txt'

    if (nrank == 0) then
     ! Open the file for writing
      if (istep==0) then
        open(unit=10, file=filename, status='old', action='write')
        write(10, '(4(ES20.10, 1X))') time, kib, eb, enst
        close(10)
      else
        open(unit=10, file=filename, status='old', position='append', action='write')
        write(10, '(4(ES20.10, 1X))') time, kib, eb, enst
        close(10)
      endif
    endif
#endif
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
  subroutine calcStats(qave,qtime,factAvg,countAvg,rho,u,v,w,ien,pre,tem,mu,ka)
    use decomp_2d
    use mod_param
    implicit none
    integer i,j,k,h,c,countAvg
    real(mytype) :: factAvg
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: rho,u,v,w,ien,pre,tem,mu,ka
    real(mytype) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, xa, za
    real(mytype) :: dilla2, dtx, dty, dtz
    real(mytype) :: cp_local, ent_local, sos_local, Pr_local
#if defined(BL) ||  defined(TGV) 
    real(mytype), dimension(:,:,:), intent(inout) :: qave,qtime
    real(mytype) :: rho_2d,u_2d,v_2d,w_2d,pre_2d,tem_2d,ien_2d,mu_2d,ka_2d,cp_2d,ent_2d,sos_2d,Pr_2d
    real(mytype) :: rhou_2d,rhov_2d,rhow_2d,rhot_2d,rhoent_2d
    real(mytype) :: uu_2d,uv_2d,uw_2d,vv_2d,vw_2d,ww_2d,rhorho_2d,temtem_2d,mumu_2d,kaka_2d,PrPr_2d,cpcp_2d
    real(mytype) :: rhouu_2d,rhouv_2d,rhouw_2d,rhovv_2d,rhovw_2d,rhoww_2d,rhouT_2d              
    real(mytype) :: tauxx_2d,tauxy_2d,tauxz_2d,tauyy_2d,tauyz_2d,tauzz_2d,qx_2d,qy_2d,qz_2d
    real(mytype) :: mudwdx_2d,mududz_2d,dwdx_2d,dudz_2d,dtdx_2d
    real(mytype), allocatable, dimension(:,:) :: tauxz_wall
#elif defined(CHA) 
    real(mytype), dimension(:,:), intent(inout) :: qave,qtime
    real(mytype) :: rho_1d,u_1d,v_1d,w_1d,pre_1d,tem_1d,ien_1d,mu_1d,ka_1d,cp_1d,ent_1d,sos_1d,Pr_1d
    real(mytype) :: rhou_1d,rhov_1d,rhow_1d,rhot_1d,rhoent_1d
    real(mytype) :: uu_1d,uv_1d,uw_1d,vv_1d,vw_1d,ww_1d,rhorho_1d,temtem_1d,mumu_1d,kaka_1d,PrPr_1d,cpcp_1d
    real(mytype) :: rhouu_1d,rhouv_1d,rhouw_1d,rhovv_1d,rhovw_1d,rhoww_1d,rhouT_1d             
    real(mytype) :: tauxx_1d,tauxy_1d,tauxz_1d,tauyy_1d,tauyz_1d,tauzz_1d,qx_1d,qy_1d,qz_1d
    real(mytype) :: mudwdx_1d,mududz_1d,dwdx_1d,dudz_1d,dtdx_1d
#endif
#if defined(BL)
    ! for boundary layer 
    !$acc parallel loop gang collapse(2) default(present)
    do k=1,xsize(3)
      do i=1,xsize(1)
        rho_2d    = 0.0_mytype
        u_2d      = 0.0_mytype
        v_2d      = 0.0_mytype
        w_2d      = 0.0_mytype
        pre_2d    = 0.0_mytype
        tem_2d    = 0.0_mytype
        ien_2d    = 0.0_mytype
        mu_2d     = 0.0_mytype
        ka_2d     = 0.0_mytype
        cp_2d     = 0.0_mytype
        ent_2d    = 0.0_mytype
        sos_2d    = 0.0_mytype
        Pr_2d     = 0.0_mytype
        rhou_2d   = 0.0_mytype
        rhov_2d   = 0.0_mytype
        rhow_2d   = 0.0_mytype
        rhot_2d   = 0.0_mytype
        rhoent_2d = 0.0_mytype
        uu_2d     = 0.0_mytype
        uv_2d     = 0.0_mytype
        uw_2d     = 0.0_mytype
        vv_2d     = 0.0_mytype
        vw_2d     = 0.0_mytype
        ww_2d     = 0.0_mytype
        rhorho_2d = 0.0_mytype
        temtem_2d = 0.0_mytype
        mumu_2d   = 0.0_mytype
        kaka_2d   = 0.0_mytype
        PrPr_2d   = 0.0_mytype
        cpcp_2d   = 0.0_mytype
        rhouu_2d  = 0.0_mytype
        rhouv_2d  = 0.0_mytype
        rhouw_2d  = 0.0_mytype
        rhovv_2d  = 0.0_mytype
        rhovw_2d  = 0.0_mytype
        rhoww_2d  = 0.0_mytype
        rhouT_2d  = 0.0_mytype
        tauxx_2d  = 0.0_mytype
        tauxy_2d  = 0.0_mytype
        tauxz_2d  = 0.0_mytype
        tauyy_2d  = 0.0_mytype
        tauyz_2d  = 0.0_mytype
        tauzz_2d  = 0.0_mytype
        qx_2d     = 0.0_mytype
        qy_2d     = 0.0_mytype
        qz_2d     = 0.0_mytype
        mudwdx_2d = 0.0_mytype
        mududz_2d = 0.0_mytype
        dwdx_2d   = 0.0_mytype
        dudz_2d   = 0.0_mytype
        dtdx_2d   = 0.0_mytype
        !$acc loop reduction(+:rho_2d,u_2d,v_2d,w_2d,pre_2d,tem_2d,ien_2d,mu_2d,ka_2d,cp_2d,ent_2d,sos_2d,Pr_2d) &
        !$acc reduction(+:rhou_2d,rhov_2d,rhow_2d,rhot_2d,rhoent_2d,uu_2d,uv_2d,uw_2d,vv_2d,vw_2d,ww_2d) &
        !$acc reduction(+:rhouu_2d,rhouv_2d,rhouw_2d,rhovv_2d,rhovw_2d,rhoww_2d,rhouT_2d) &
        !$acc reduction(+:rhorho_2d,temtem_2d,mumu_2d,kaka_2d,PrPr_2d,cpcp_2d) &
        !$acc reduction(+:tauxx_2d,tauxy_2d,tauxz_2d,tauyy_2d,tauyz_2d,tauzz_2d,qx_2d,qy_2d,qz_2d) &
        !$acc reduction(+:dwdx_2d,dudz_2d,dtdx_2d,mudwdx_2d,mududz_2d) &
        !$acc private(dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz,dtx,dty,dtz,cp_local,ent_local,sos_local,Pr_local) 
        do j=1,xsize(2)
          call calcCpH(rho(i,j,k),ien(i,j,k),cp_local,ent_local)
          call calcSOS(rho(i,j,k),ien(i,j,k),sos_local)
          rho_2d    = rho_2d + rho(i,j,k)
          u_2d      = u_2d + u(i,j,k)
          v_2d      = v_2d + v(i,j,k)
          w_2d      = w_2d + w(i,j,k)
          pre_2d    = pre_2d + pre(i,j,k)
          tem_2d    = tem_2d + tem(i,j,k)
          ien_2d    = ien_2d + ien(i,j,k)
          mu_2d     = mu_2d + mu(i,j,k)
          ka_2d     = ka_2d + ka(i,j,k)
          cp_2d     = cp_2d + cp_local
          Pr_local  = cp_local*mu(i,j,k)/ka(i,j,k)
          Pr_2d     = Pr_2d  + Pr_local
          ent_2d    = ent_2d + ent_local
          sos_2d    = sos_2d + sos_local
          rhou_2d   = rhou_2d + rho(i,j,k)*u(i,j,k)
          rhov_2d   = rhov_2d + rho(i,j,k)*v(i,j,k)
          rhow_2d   = rhow_2d + rho(i,j,k)*w(i,j,k)
          rhot_2d   = rhot_2d + rho(i,j,k)*tem(i,j,k)
          rhoent_2d = rhoent_2d + rho(i,j,k)*ent_local
          uu_2d     = uu_2d + u(i,j,k)*u(i,j,k)
          uv_2d     = uv_2d + u(i,j,k)*v(i,j,k)
          uw_2d     = uw_2d + u(i,j,k)*w(i,j,k)
          vv_2d     = vv_2d + v(i,j,k)*v(i,j,k)
          vw_2d     = vw_2d + v(i,j,k)*w(i,j,k)
          ww_2d     = ww_2d + w(i,j,k)*w(i,j,k)
          rhorho_2d = rhorho_2d + rho(i,j,k)*rho(i,j,k)
          temtem_2d = temtem_2d + tem(i,j,k)*tem(i,j,k)
          mumu_2d   = mumu_2d + mu(i,j,k)*mu(i,j,k)
          kaka_2d   = kaka_2d + ka(i,j,k)*ka(i,j,k)
          PrPr_2d   = PrPr_2d + Pr_local*Pr_local
          cpcp_2d   = cpcp_2d + cp_local*cp_local
          rhouu_2d  = rhouu_2d + rho(i,j,k)*u(i,j,k)*u(i,j,k)
          rhouv_2d  = rhouv_2d + rho(i,j,k)*u(i,j,k)*v(i,j,k)
          rhouw_2d  = rhouw_2d + rho(i,j,k)*u(i,j,k)*w(i,j,k)
          rhovv_2d  = rhovv_2d + rho(i,j,k)*v(i,j,k)*v(i,j,k)
          rhovw_2d  = rhovw_2d + rho(i,j,k)*v(i,j,k)*w(i,j,k)
          rhoww_2d  = rhoww_2d + rho(i,j,k)*w(i,j,k)*w(i,j,k)
          rhouT_2d  = rhouT_2d + rho(i,j,k)*u(i,j,k)*tem(i,j,k)
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
          dilla2     = dux + dvy + dwz
          tauxx_2d   = tauxx_2d + mu(i,j,k)*(2.0_mytype*dux - 2.0_mytype/3.0_mytype*dilla2) 
          tauxy_2d   = tauxy_2d + mu(i,j,k)*(duy + dvx)
          tauxz_2d   = tauxz_2d + mu(i,j,k)*(duz + dwx)
          tauyy_2d   = tauyy_2d + mu(i,j,k)*(2.0_mytype*dvy - 2.0_mytype/3.0_mytype*dilla2)
          tauyz_2d   = tauyz_2d + mu(i,j,k)*(dvz + dwy)
          tauzz_2d   = tauzz_2d + mu(i,j,k)*(2.0_mytype*dwz - 2.0_mytype/3.0_mytype*dilla2)
          qx_2d      = qx_2d    + ka(i,j,k)*dtx
          qy_2d      = qy_2d    + ka(i,j,k)*dty
          qz_2d      = qz_2d    + ka(i,j,k)*dtz
          mudwdx_2d  = mudwdx_2d  + mu(i,j,k)*dwx
          mududz_2d  = mududz_2d  + mu(i,j,k)*duz
          dwdx_2d    = dwdx_2d  + dwx
          dudz_2d    = dudz_2d  + duz   
          dtdx_2d    = dtdx_2d  + dtx 
        enddo
      ! primary variables
      qave(i,k,1)  = qave(i,k,1) + rho_2d  
      qave(i,k,2)  = qave(i,k,2) + u_2d  
      qave(i,k,3)  = qave(i,k,3) + v_2d 
      qave(i,k,4)  = qave(i,k,4) + w_2d
      qave(i,k,5)  = qave(i,k,5) + pre_2d
      qave(i,k,6)  = qave(i,k,6) + tem_2d
      qave(i,k,7)  = qave(i,k,7) + ien_2d
      qave(i,k,8)  = qave(i,k,8) + mu_2d
      qave(i,k,9)  = qave(i,k,9) + ka_2d
      qave(i,k,10) = qave(i,k,10) + cp_2d
      qave(i,k,11) = qave(i,k,11) + ent_2d
      qave(i,k,12) = qave(i,k,12) + sos_2d
      qave(i,k,13) = qave(i,k,13) + Pr_2d
      ! Favre products
      qave(i,k,14) = qave(i,k,14) + rhou_2d
      qave(i,k,15) = qave(i,k,15) + rhov_2d
      qave(i,k,16) = qave(i,k,16) + rhow_2d
      qave(i,k,17) = qave(i,k,17) + rhot_2d
      qave(i,k,18) = qave(i,k,18) + rhoent_2d
      ! double products
      qave(i,k,19) = qave(i,k,19) + uu_2d
      qave(i,k,20) = qave(i,k,20) + uv_2d
      qave(i,k,21) = qave(i,k,21) + uw_2d
      qave(i,k,22) = qave(i,k,22) + vv_2d
      qave(i,k,23) = qave(i,k,23) + vw_2d
      qave(i,k,24) = qave(i,k,24) + ww_2d
      qave(i,k,25) = qave(i,k,25) + rhorho_2d
      qave(i,k,26) = qave(i,k,26) + temtem_2d
      qave(i,k,27) = qave(i,k,27) + mumu_2d
      qave(i,k,28) = qave(i,k,28) + kaka_2d
      qave(i,k,29) = qave(i,k,29) + PrPr_2d
      qave(i,k,30) = qave(i,k,30) + cpcp_2d
      ! Favre double products
      qave(i,k,31) = qave(i,k,31) + rhouu_2d 
      qave(i,k,32) = qave(i,k,32) + rhouv_2d 
      qave(i,k,33) = qave(i,k,33) + rhouw_2d 
      qave(i,k,34) = qave(i,k,34) + rhovv_2d 
      qave(i,k,35) = qave(i,k,35) + rhovw_2d 
      qave(i,k,36) = qave(i,k,36) + rhoww_2d 
      qave(i,k,37) = qave(i,k,37) + rhouT_2d
      ! stress tensor
      qave(i,k,38) = qave(i,k,38) + tauxx_2d
      qave(i,k,39) = qave(i,k,39) + tauxy_2d 
      qave(i,k,40) = qave(i,k,40) + tauxz_2d 
      qave(i,k,41) = qave(i,k,41) + tauyy_2d 
      qave(i,k,42) = qave(i,k,42) + tauyz_2d 
      qave(i,k,43) = qave(i,k,43) + tauzz_2d
      ! heat flux 
      qave(i,k,44) = qave(i,k,44) + qx_2d
      qave(i,k,45) = qave(i,k,45) + qy_2d
      qave(i,k,46) = qave(i,k,46) + qz_2d
      ! extras
      qave(i,k,47) = qave(i,k,47) + dwdx_2d
      qave(i,k,48) = qave(i,k,48) + dudz_2d
      qave(i,k,49) = qave(i,k,49) + dtdx_2d
      qave(i,k,50) = qave(i,k,50) + mudwdx_2d
      qave(i,k,51) = qave(i,k,51) + mududz_2d
    enddo
    enddo

    ! for boundary layer
    !$acc parallel loop gang collapse(2) private(duz,dwx) 
    do k=1,xsize(3)
      do j=1,xsize(2)
        duz = 0.0_mytype
        dwx = 0.0_mytype
        xa = xp(1)
        za = zp(k)
        !$acc loop seq
        do c = 1, nStencilVisc
          duz = duz + visc_ddz(c)*(u(1,j,k+c) - u(1,j,k-c))*za
          dwx = dwx + visc_ddx(c)*(w(1+c,j,k) - w(1-c,j,k))*xa
        enddo
        qtime(j,k,1)  = qtime(j,k,1) + mu(1,j,k)*(duz + dwx)
      enddo
    enddo
    factAvg = factAvg + 1.0_mytype*xsize(2)

#elif defined(CHA)
    ! for channel
    !$acc parallel loop gang default(present)
    do i=1,xsize(1)
      rho_1d    = 0.0_mytype
      u_1d      = 0.0_mytype
      v_1d      = 0.0_mytype
      w_1d      = 0.0_mytype
      pre_1d    = 0.0_mytype
      tem_1d    = 0.0_mytype
      ien_1d    = 0.0_mytype
      mu_1d     = 0.0_mytype
      ka_1d     = 0.0_mytype
      cp_1d     = 0.0_mytype
      ent_1d    = 0.0_mytype
      sos_1d    = 0.0_mytype
      Pr_1d     = 0.0_mytype
      rhou_1d   = 0.0_mytype
      rhov_1d   = 0.0_mytype
      rhow_1d   = 0.0_mytype
      rhot_1d   = 0.0_mytype
      rhoent_1d = 0.0_mytype
      uu_1d     = 0.0_mytype
      uv_1d     = 0.0_mytype
      uw_1d     = 0.0_mytype
      vv_1d     = 0.0_mytype
      vw_1d     = 0.0_mytype
      ww_1d     = 0.0_mytype
      rhorho_1d = 0.0_mytype
      temtem_1d = 0.0_mytype
      mumu_1d   = 0.0_mytype
      kaka_1d   = 0.0_mytype
      PrPr_1d   = 0.0_mytype
      cpcp_1d   = 0.0_mytype
      rhouu_1d  = 0.0_mytype
      rhouv_1d  = 0.0_mytype
      rhouw_1d  = 0.0_mytype
      rhovv_1d  = 0.0_mytype
      rhovw_1d  = 0.0_mytype
      rhoww_1d  = 0.0_mytype
      rhouT_1d  = 0.0_mytype
      tauxx_1d  = 0.0_mytype
      tauxy_1d  = 0.0_mytype
      tauxz_1d  = 0.0_mytype
      tauyy_1d  = 0.0_mytype
      tauyz_1d  = 0.0_mytype
      tauzz_1d  = 0.0_mytype
      qx_1d     = 0.0_mytype
      qy_1d     = 0.0_mytype
      qz_1d     = 0.0_mytype
      mudwdx_1d = 0.0_mytype
      mududz_1d = 0.0_mytype
      dwdx_1d   = 0.0_mytype
      dudz_1d   = 0.0_mytype
      dtdx_1d   = 0.0_mytype
      !$acc loop collapse(2) reduction(+:rho_1d,u_1d,v_1d,w_1d,pre_1d,tem_1d,ien_1d,mu_1d,ka_1d,cp_1d,ent_1d,sos_1d,Pr_1d) &
      !$acc reduction(+:rhou_1d,rhov_1d,rhow_1d,rhot_1d,rhoent_1d,uu_1d,uv_1d,uw_1d,vv_1d,vw_1d,ww_1d) &
      !$acc reduction(+:rhorho_1d,temtem_1d,mumu_1d,kaka_1d,PrPr_1d,cpcp_1d) &
      !$acc reduction(+:rhouu_1d,rhouv_1d,rhouw_1d,rhovv_1d,rhovw_1d,rhoww_1d,rhouT_1d) &
      !$acc reduction(+:tauxx_1d,tauxy_1d,tauxz_1d,tauyy_1d,tauyz_1d,tauzz_1d,qx_1d,qy_1d,qz_1d) &
      !$acc reduction(+:dwdx_1d,dudz_1d,dtdx_1d,mudwdx_1d,mududz_1d) &
      !$acc private(dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz,dtx,dty,dtz,cp_local,ent_local,sos_local,Pr_local) 
      do k=1,xsize(3)
        do j=1,xsize(2)
          call calcCpH(rho(i,j,k),ien(i,j,k),cp_local,ent_local)
          call calcSOS(rho(i,j,k),ien(i,j,k),sos_local)
          rho_1d   = rho_1d + rho(i,j,k)
          u_1d     = u_1d + u(i,j,k)
          v_1d     = v_1d + v(i,j,k)
          w_1d     = w_1d + w(i,j,k)
          pre_1d   = pre_1d + pre(i,j,k)
          tem_1d   = tem_1d + tem(i,j,k)
          ien_1d   = ien_1d + ien(i,j,k)
          mu_1d    = mu_1d + mu(i,j,k)
          ka_1d    = ka_1d + ka(i,j,k)
          cp_1d    = cp_1d + cp_local
          Pr_local = cp_local*mu(i,j,k)/ka(i,j,k)
          Pr_1d    = Pr_1d  + Pr_local
          ent_1d   = ent_1d + ent_local
          sos_1d   = sos_1d + sos_local
          rhou_1d  = rhou_1d + rho(i,j,k)*u(i,j,k)
          rhov_1d  = rhov_1d + rho(i,j,k)*v(i,j,k)
          rhow_1d  = rhow_1d + rho(i,j,k)*w(i,j,k)
          rhot_1d  = rhot_1d + rho(i,j,k)*tem(i,j,k)
          uu_1d    = uu_1d + u(i,j,k)*u(i,j,k)
          uv_1d    = uv_1d + u(i,j,k)*v(i,j,k)
          uw_1d    = uw_1d + u(i,j,k)*w(i,j,k)
          vv_1d    = vv_1d + v(i,j,k)*v(i,j,k)
          vw_1d    = vw_1d + v(i,j,k)*w(i,j,k)
          ww_1d    = ww_1d + w(i,j,k)*w(i,j,k)
          rhorho_1d = rhorho_1d + rho(i,j,k)*rho(i,j,k)
          temtem_1d = temtem_1d + tem(i,j,k)*tem(i,j,k)
          mumu_1d   = mumu_1d + mu(i,j,k)*mu(i,j,k)
          kaka_1d   = kaka_1d + ka(i,j,k)*ka(i,j,k)
          PrPr_1d   = temtem_1d + Pr_local*Pr_local
          cpcp_1d   = temtem_1d + cp_local*cp_local
          rhouu_1d = rhouu_1d + rho(i,j,k)*u(i,j,k)*u(i,j,k)
          rhouv_1d = rhouv_1d + rho(i,j,k)*u(i,j,k)*v(i,j,k)
          rhouw_1d = rhouw_1d + rho(i,j,k)*u(i,j,k)*w(i,j,k)
          rhovv_1d = rhovv_1d + rho(i,j,k)*v(i,j,k)*v(i,j,k)
          rhovw_1d = rhovw_1d + rho(i,j,k)*v(i,j,k)*w(i,j,k)
          rhoww_1d = rhoww_1d + rho(i,j,k)*w(i,j,k)*w(i,j,k)
          rhouT_1d = rhouT_1d + rho(i,j,k)*u(i,j,k)*tem(i,j,k)
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
          dilla2     = dux + dvy + dwz
          tauxx_1d   = tauxx_1d + mu(i,j,k)*(2.0_mytype*dux - 2.0_mytype/3.0_mytype*dilla2) 
          tauxy_1d   = tauxy_1d + mu(i,j,k)*(duy + dvx)
          tauxz_1d   = tauxz_1d + mu(i,j,k)*(duz + dwx)
          tauyy_1d   = tauyy_1d + mu(i,j,k)*(2.0_mytype*dvy - 2.0_mytype/3.0_mytype*dilla2)
          tauyz_1d   = tauyz_1d + mu(i,j,k)*(dvz + dwy)
          tauzz_1d   = tauzz_1d + mu(i,j,k)*(2.0_mytype*dwz - 2.0_mytype/3.0_mytype*dilla2)
          qx_1d      = qx_1d    + ka(i,j,k)*dtx
          qy_1d      = qy_1d    + ka(i,j,k)*dty
          qz_1d      = qz_1d    + ka(i,j,k)*dtz
          mudwdx_1d  = mudwdx_1d  + mu(i,j,k)*dwx
          mududz_1d  = mududz_1d  + mu(i,j,k)*duz
          dwdx_1d    = dwdx_1d  + dwx  
          dudz_1d    = dudz_1d  + duz 
          dtdx_1d    = dtdx_1d  + dtx 
        enddo
      enddo
      ! primary variables
      qave(i,1)  = qave(i,1) + rho_1d  
      qave(i,2)  = qave(i,2) + u_1d  
      qave(i,3)  = qave(i,3) + v_1d 
      qave(i,4)  = qave(i,4) + w_1d
      qave(i,5)  = qave(i,5) + pre_1d
      qave(i,6)  = qave(i,6) + tem_1d
      qave(i,7)  = qave(i,7) + ien_1d
      qave(i,8)  = qave(i,8) + mu_1d
      qave(i,9)  = qave(i,9) + ka_1d
      qave(i,10) = qave(i,10) + cp_1d
      qave(i,11) = qave(i,11) + ent_1d
      qave(i,12) = qave(i,12) + sos_1d
      qave(i,13) = qave(i,13) + Pr_1d
      ! favre products
      qave(i,14) = qave(i,14) + rhou_1d
      qave(i,15) = qave(i,15) + rhov_1d
      qave(i,16) = qave(i,16) + rhow_1d
      qave(i,17) = qave(i,17) + rhot_1d
      qave(i,18) = qave(i,18) + rhoent_1d
      ! double products
      qave(i,19) = qave(i,19) + uu_1d
      qave(i,20) = qave(i,20) + uv_1d
      qave(i,21) = qave(i,21) + uw_1d
      qave(i,22) = qave(i,22) + vv_1d
      qave(i,23) = qave(i,23) + vw_1d
      qave(i,24) = qave(i,24) + ww_1d
      qave(i,25) = qave(i,25) + rhorho_1d
      qave(i,26) = qave(i,26) + temtem_1d
      qave(i,27) = qave(i,27) + mumu_1d
      qave(i,28) = qave(i,28) + kaka_1d
      qave(i,29) = qave(i,29) + PrPr_1d
      qave(i,30) = qave(i,30) + cpcp_1d
      ! favre double products
      qave(i,31) = qave(i,31) + rhouu_1d
      qave(i,32) = qave(i,32) + rhouv_1d
      qave(i,33) = qave(i,33) + rhouw_1d
      qave(i,34) = qave(i,34) + rhovv_1d
      qave(i,35) = qave(i,35) + rhovw_1d
      qave(i,36) = qave(i,36) + rhoww_1d
      qave(i,37) = qave(i,37) + rhouT_1d
      ! stress tensor
      qave(i,38) = qave(i,38) + tauxx_1d
      qave(i,39) = qave(i,39) + tauxy_1d 
      qave(i,40) = qave(i,40) + tauxz_1d 
      qave(i,41) = qave(i,41) + tauyy_1d 
      qave(i,42) = qave(i,42) + tauyz_1d 
      qave(i,43) = qave(i,43) + tauzz_1d
      ! heat flux 
      qave(i,44) = qave(i,44) + qx_1d
      qave(i,45) = qave(i,45) + qy_1d
      qave(i,46) = qave(i,46) + qz_1d
      ! extras
      qave(i,47) = qave(i,47) + dwdx_1d
      qave(i,48) = qave(i,48) + dudz_1d
      qave(i,49) = qave(i,49) + dtdx_1d
      qave(i,50) = qave(i,50) + mudwdx_1d
      qave(i,51) = qave(i,51) + mududz_1d
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
  subroutine loadRestart(istep,time,dp,rho,u,v,w,ien,nHaloIn,part)
    use mpi
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,ierr,nHaloIn
    character*7 cha
    real(mytype), dimension(1-nHaloIn:,1-nHaloIn:,1-nHaloIn:) :: rho,u,v,w,ien
    real(mytype), allocatable, dimension(:,:,:) :: tmp
    real(mytype), dimension(6) :: infoRestart
    real(mytype) :: time, dp
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
    call decomp_2d_read_scalar(fh,disp,6,infoRestart)
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
    dp = infoRestart(5)
    time = infoRestart(6)
    if(nrank==0) write(stdout,'(A,E12.5,A,I0)') 'The simulation restarts at t: ',time, ' and at step: ', istep
    if(nrank==0) write(stdout,*) ''
  end subroutine


! I/O save restart
  subroutine saveRestart(istep,time,dp,rho,u,v,w,ien,nHaloIn,part,interpol)
    use mpi
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,ierr,nHaloIn
    character*7 cha
    real(mytype), dimension(1-nHaloIn:,1-nHaloIn:,1-nHaloIn:) :: rho,u,v,w,ien
    real(mytype), allocatable, dimension(:,:,:) :: tmp
    real(mytype), dimension(6) :: infoRestart
    real(mytype) :: time, dp
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
    infoRestart(5) = dp
    infoRestart(6) = time 
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
    call decomp_2d_write_scalar(fh,disp,6,infoRestart)
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
  subroutine saveStats(part,istep,dt,qave,qtime,factAvg,countAvg)
    use decomp_2d
    use decomp_2d_io
    use mod_param
    implicit none
    integer :: istep,countAvg
    character*7 cha
    real(mytype) :: dt, factAvg, factAvg_inv, countAvg_inv
    TYPE (DECOMP_INFO) :: part
    logical :: exist
    real(mytype), allocatable, dimension(:,:,:) :: tmpPlane
#if defined(BL) || defined(TGV)
    ! for boundary layer
    real(mytype), dimension(:,:,:), intent(in) :: qave
    real(mytype), dimension(:,:,:), intent(in) :: qtime
#elif defined(CHA)
    real(mytype), dimension(:,:), intent(in) :: qave
    real(mytype), dimension(:,:), intent(in) :: qtime
#endif
    allocate(tmpPlane(part%xsz(1),part%xsz(2),part%xsz(3)))
    tmpPlane = 0.0_mytype
    write(cha,'(I0.7)') istep
    factAvg_inv = 1.0_mytype/factAvg
    countAvg_inv = 1.0_mytype/countAvg
    ! write stats planes in /output/stats/ 
    ! for boundary layer
#if defined(BL)
    call write_all_qave_51()
    ! tau_xz_wall (1-time)
    tmpPlane(1,:,:) = qtime(:,:,1) * countAvg_inv
    call decomp_2d_write_plane(1, tmpPlane, 1, 1, '.', &
         'output/stats/tauxz_wall_time.'//cha//'.bin', 'dummy')
#elif defined(CHA)
    call write_all_qave_51()
#endif
    if (nrank == 0) then
#if defined(BL)
      write(stdout,'(A, I0)') 'Boundary layer: writing time- and span-averaged statistics at final step: ', istep
#elif defined(CHA)
      write(stdout,'(A, I0)') 'Channel: writing time-, span- and stream-averaged statistics at final step: ', istep
#endif
      inquire(file="output/stats/stats_info.txt", exist=exist)
      if (exist) then
        open(11,file='output/stats/stats_info.txt',status="old",position="append",action="write")
      else
        open(11,file='output/stats/stats_info.txt',status="new",action="write")
      endif
      write(11,'(A, I0)')   'step: ', istep
      write(11,'(A, I0)')   'number of averaged files: ', countAvg
      write(11,'(A, I0)')   'average factor: ', int(factAvg)
      write(11,'(A, E17.6)') 'time from begin averaging: ', dt*(istep - saveStatsAfter)
      close(11)
    endif
    deallocate(tmpPlane)

  contains

    subroutine write_stat_plane(idx, base)
      integer, intent(in) :: idx
      character(len=*), intent(in) :: base
      character(len=256) :: fname
      integer :: u
      fname = 'output/stats/'//trim(base)//cha//'.bin'
#if defined(BL)
      tmpPlane(:,1,:) = qave(:,:,idx) * factAvg_inv
      call decomp_2d_write_plane(1, tmpPlane, 2, 1, '.', trim(fname), 'dummy')
#elif defined(CHA)
      if (nrank == 0) then
        u = 99
        open(u, file=trim(fname), form='unformatted', access='stream', &
             status='replace', action='write')
        write(u) qave(:,idx) * factAvg_inv
        close(u)
      endif
#endif 
    end subroutine write_stat_plane

    subroutine write_all_qave_51()
      integer, parameter :: nvar = 51
      integer :: n
      character(len=32), parameter :: tag(nvar) = [ character(len=32) :: &
        'r_avg.',      'u_avg.',      'v_avg.',      'w_avg.',      &
        'p_avg.',      't_avg.',      'ien_avg.',    'mu_avg.',     &
        'ka_avg.',     'cp_avg.',     'ent_avg.',    'sos_avg.',    &
        'Pr_avg.',     'ru_avg.',     'rv_avg.',     'rw_avg.',     &
        'rt_avg.',     'rent_avg.',   'uu_avg.',     'uv_avg.',     &
        'uw_avg.',     'vv_avg.',     'vw_avg.',     'ww_avg.',     &
        'rhorho_avg.', 'temtem_avg.', 'mumu_avg.',   'kaka_avg.',   &
        'PrPr_avg.',   'cpcp_avg.',   'ruu_avg.',    'ruv_avg.',    &
        'ruw_avg.',    'rvv_avg.',    'rvw_avg.',    'rww_avg.',    &
        'rutem_avg.',  'tauxx_avg.',  'tauxy_avg.',  'tauxz_avg.',  &
        'tauyy_avg.',  'tauyz_avg.',  'tauzz_avg.',  'qx_avg.',     &
        'qy_avg.',     'qz_avg.',     'dwdx_avg.',   'dudz_avg.',   &
        'dtdx_avg.',   'mudwdx_avg.', 'mududz_avg.' ]
      do n = 1, nvar
        call write_stat_plane(n, tag(n))
      end do
    end subroutine write_all_qave_51
  end subroutine saveStats

end module
