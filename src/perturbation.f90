! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_perturbation
  use mod_param
  use mod_grid
  use mod_math
  implicit none
  integer :: kC, LP, kSt, kEn, n3dmode
  real(mytype), allocatable, dimension(:,:) :: dudt_pert
  real(mytype), allocatable, dimension(:) :: modes3d
  real(mytype), dimension(1:2) :: pert_omeg
  real(mytype) :: lamda, pert_zStart, pert_zEnd, freq1, freq2, amp1, amp2
  real(mytype) :: pert_zMid_new, pert_zStart_new, pert_zEnd_new, pert_zLen_new, pert_Restart_new, pert_Remid_new, pert_Reend_new
contains
  ! initialization perturbation module
  subroutine init_pert()
    use decomp_2d
    use mod_param
    use mod_grid
    implicit none
    integer :: ierr,j,m 
    real(mytype) :: y_glob
    ! dudt_pert is the time variation of the wall-normal velocity
    allocate(dudt_pert(1-nHalo:xsize(2)+nHalo, 1-nHalo:xsize(3)+nHalo))
    allocate(modes3d(xsize(2)))
    dudt_pert = 0.0_mytype
    !$acc enter data copyin(dudt_pert)
    if (pert_calc == 1) then
      ! placement of the disturbance strip
      pert_zMid = pert_ReMid**2/Re-zStartDNS
      pert_zStart = pert_zMid-pert_zLen/2
      pert_zEnd = pert_zMid+pert_zLen/2
      ! find index for starting, mid, and ending point
      kC  = minloc(abs(pert_zMid-z_global),1)
      kSt = minloc(abs(pert_zStart-z_global),1)
      kEn = minloc(abs(pert_zEnd-z_global),1)
      ! recalculate the perturbation location
      pert_zMid_new = z_global(kC)
      pert_zStart_new = z_global(kSt)
      pert_zEnd_new = z_global(kEn)
      pert_zLen_new = pert_zEnd_new - pert_zStart_new
      pert_Restart_new = (Re)**0.5_mytype*( (pert_zStart_new+zStartDNS)**0.5_mytype )
      pert_Remid_new = (Re)**0.5_mytype*( (pert_zMid_new+zStartDNS)**0.5_mytype )
      pert_Reend_new = (Re)**0.5_mytype*( (pert_zEnd_new+zStartDNS)**0.5_mytype )
      ! extraction of the frequencies from config.h
      freq1=pert_F(1)
      freq2=pert_F(2)
      omega2 = Re*freq2
      pert_omeg=(/omega1,omega2/)
      ! extraction of the amplitudes from config.h
      amp1=pert_ampl(1)
      amp2=pert_ampl(2)
      ! if the perturbation is placed inside the sponge, then stop the simulation
      if (pert_Restart_new .lt. spInLen_Reend) then
        if (nrank == 0) write(stdout,*) "suction/blowing is inside the inlet sponge! RESET!"
          call decomp_2d_finalize
          call mpi_finalize(ierr) 
          stop 
      endif
      ! spanwise part of the disturbance
      n3dmode=size(pert_beta)
      modes3d=0.0_mytype
      do j = 1,xsize(2)
        y_glob = (j + xstart(2) - 1)*dy
        do m=1,n3dmode
          modes3d(j) = modes3d(j) + cos(pert_beta(m)*y_glob)
        enddo
      enddo
      !$acc enter data copyin(modes3d)
      ! Sayadi et al., JFM 724, 2013
      !$acc enter data copyin(kC,kEn,kSt,amp1,amp2,omega1,omega2,pert_omeg,pert_ampl,pert_beta,n3dmode)
      ! Franko & Lele, JFM 730, 2013
      !!!$acc enter data copyin(pert_omeg,pert_ampl,pert_beta,n3dmode,pert_ySig,pert_zMid_new,pert_zStart_new,pert_zEnd_new)
    endif
  end subroutine
! definition of the disturbance strip according to Sayadi et al., JFM 724, 2013
subroutine perturbationSayadi(u,time, i_index, dudt_pert)
  use decomp_2d
  use mod_param
  implicit none
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: u
  real(mytype), dimension(1-nHalo:,1-nHalo:) :: dudt_pert
  real(mytype) :: time, ksi, g, kappa, alpha, beta, y_glob, fact_u_2D, fact_u_3D, fact_du_2D, fact_du_3D
  integer :: i,j,k,c, kk, i_index, m
  !$acc parallel loop default(present) private(kk,alpha,beta,ksi,g,y_glob,fact_u_2D,fact_u_3D,fact_du_2D,fact_du_3D) async(1)
  do k = kSt,kEn
    if ((k >= xstart(3)) .and. (k <= xend(3))) then
      kk = (1.0_mytype*k-xstart(3)+1)
      if (k <= kC) then
        kappa = 1.0_mytype
        alpha = k - kSt
        beta = kC - kSt
      else
        kappa = -1.0_mytype
        alpha = kEn - k
        beta = kEn - kC
      endif
      ksi = alpha/beta
      g = (15.1875_mytype*ksi**5) - (35.4375_mytype*ksi**4) + (20.25_mytype*ksi**3)
      fact_u_2D = amp1*kappa*g*sin(omega1*time)
      fact_u_3D = amp2*kappa*g*sin(omega2*time)
      fact_du_2D = amp1*kappa*g*omega1*cos(omega1*time)
      fact_du_3D = amp2*kappa*g*omega2*cos(omega2*time)
      !$acc loop seq
      do j = 1,xsize(2)
        u(i_index,j,kk) = fact_u_3D*modes3d(j)+fact_u_2D
        dudt_pert(j,kk) = fact_du_3D*modes3d(j)+fact_du_2D  
      enddo
    endif
  enddo
end subroutine
! definition of the disturbance strip according to Franko & Lele, JFM 730, 2013
subroutine perturbationFrankoLele(u,time, i_index,dudt_pert)
  use decomp_2d
  use mod_param
  implicit none
  real(mytype), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: u
  real(mytype), dimension(1-nHalo:,1-nHalo:) :: dudt_pert
  real(mytype), dimension(1:n3dmode) :: sum_u, sum_du
  real(mytype) :: time, zz, fz, yy, fy
  integer :: i,j,k,m,i_index
  !$acc parallel loop default(present) private(yy,zz,fy,fz,n3dmode,sum_u,sum_du) async(1)
  do k=1,xsize(3)
    zz = (k + xstart(3) - 2)*dz
    if ((zz > pert_zStart_new) .and. (zz < pert_zEnd_new)) then
      fz = exp(-(zz-pert_zMid_new)**2/(2.0_mytype*pert_zSig**2))
      !$acc loop seq
      do j = 1,xsize(2)
        yy = (j + xstart(2) - 1)*dy
        fy = 1.0_mytype + 0.1_mytype*(  exp(-((yy-len_y/2.0_mytype-pert_ySig)/pert_ySig)**2) &
                        - exp(-((yy-len_y/2.0_mytype+pert_ySig)/pert_ySig)**2))
        !$acc loop seq
        do m=1,n3dmode
          sum_u(m) = pert_ampl(m)*sin(pert_omeg(m)*time - pert_beta(m)*yy)
          sum_du(m) = pert_omeg(m)*pert_ampl(m)*cos(pert_omeg(m)*time - pert_beta(m)*yy)
        enddo 
      u(i_index,j,k) = sum(sum_u)*fz*fy
      dudt_pert(j,k) = sum(sum_du)*fz*fy
      enddo
    endif
  enddo 
end subroutine
! print perturbation settings
subroutine print_pertBC()
  use decomp_2d
  use mod_param
  use mod_eos
  use mod_halo
  use mod_finitediff
  implicit none
  if (nrank == 0) then
    write(stdout,* ) 'Perturbation: blowing/suction'
    write(stdout,'(A)') 'o--------------------------------------------------o'
    write(stdout,'(A, I10)') 'pert_calc:                            ', pert_calc
    if (pert_calc == 0) then
      write(stdout,'(A)') 'blowing/suction off!                                '
    elseif (pert_calc == 1) then
      write(stdout,'(A)') 'pert is outside sponge, go ahead!                 '
      write(stdout,'(A, F10.4)') 'new pert_ReMid:                       ',pert_Remid_new
      write(stdout,'(A, F10.4)') 'new pert_zMid =                       ',pert_zMid_new
      write(stdout,'(A, F10.4)') 'new pert_zLen =                       ',pert_zLen_new
      write(stdout,'(A, F10.4)') 'new pert_ReStart =                    ',pert_Restart_new
      write(stdout,'(A, F10.4)') 'new pert_ReEnd =                      ',pert_Reend_new
      write(stdout,* ) 
      if (amp2 == 0.0_mytype) then
        write(stdout,* ) 'Only 2D waves!'
        write(stdout,'(A, F10.4)') '2D dimensionless frequency (x 10^6):  ',freq1*1e6_mytype
        write(stdout,'(A, F10.4)') '2D angular frequency:                 ',omega1
        write(stdout,'(A, F10.4)') '2D wave amplitude:                    ',amp1
      else
        write(stdout,* ) '2D + oblique waves!'
        write(stdout,'(A, F10.4)') '2D dimensionless frequency (x 10^6):  ',freq1*1e6_mytype
        write(stdout,'(A, F10.4)') '2D angular frequency:                 ',omega1
        write(stdout,'(A, F10.4)') '2D wave amplitude:                    ',amp1
        write(stdout,'(A, F10.4)') '3D dimensionless frequency (x 10^6):  ',freq2*1e6_mytype
        write(stdout,'(A, F10.4)') '3D angular frequency:                 ',omega2
        write(stdout,'(A, F10.4)') '3D wave amplitude:                    ',amp2
        write(stdout,'(A, I10)') 'Number of spanwise modes:             ',n3dmode
        write(stdout,'(A, F10.4)') 'Fundamental spanwise wavenumber:      ',beta0
        if (n3dmode .gt. 1) then
          write(stdout,'(A, F10.4)') 'Spanwise wavenumbers:                 ',pert_beta
        endif
        write(stdout,* ) 
      endif
    endif
    write(stdout,'(A)') 'o--------------------------------------------------o'
    write(stdout,* ) 
  endif
  end subroutine
end module mod_perturbation
