! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_halo
  use MPI
  use mod_param, only:nHalo,BC_inl_rescale
  use decomp_2d 
  use mod_timer
  implicit none
  type rk
    integer :: myid,kp,km,jp,jm
    logical :: inlet,outlet
    integer :: kp_rcy,km_rcy
    logical :: recycle
  end type rk
  type(rk) :: neigh
  integer :: krcy
  character(len=30) :: char_buff_j
  ! definition of the transfer data arrays depending on the architecture
#if defined(_OPENACC) && !defined(_GPU_DIRECT)
  real(mytype), allocatable, dimension(:,:,:,:) :: buff_send9_j, buff_send9_k, buff_recv9_j, buff_recv9_k
  real(mytype), allocatable, dimension(:,:,:) :: buff_send1_j, buff_send1_k, buff_recv1_j, buff_recv1_k
  real(mytype), allocatable, dimension(:,:,:,:) :: send_rcy, recv_rcy
#elif defined(_OPENACC)
  real(mytype), allocatable, dimension(:,:,:,:) :: buff_send9_j, buff_send9_k, buff_recv9_j, buff_recv9_k
  real(mytype), allocatable, dimension(:,:,:) :: buff_send1_j, buff_send1_k, buff_recv1_j, buff_recv1_k
  real(mytype), allocatable, dimension(:,:,:,:) :: send_rcy, recv_rcy
#elif !defined(_OPENACC)
  real(mytype), allocatable, dimension(:,:,:) :: senjP,rcvjP,senjM,rcvjM
  real(mytype), allocatable, dimension(:,:,:) :: senkP,rcvkP,senkM,rcvkM
  real(mytype), allocatable, dimension(:,:,:,:) :: send_rcy, recv_rcy
#endif
contains
! initialization of the MPI routine
  subroutine comm_init(myid,p_row,p_col,comm_cart,xs,imax,jmax,kmax)
    implicit none
    integer               :: p_row,p_col,imax,jmax,kmax
    integer               :: rank_kp,rank_km,rank_jp,rank_jm,comm_cart,ierr,myid
    integer, dimension(2) :: dimens,coordr,coordn,coords,coordt,coordb
    logical, dimension(2) :: period
    integer, dimension(3) :: xs
#if defined(_OPENACC)
    ! sending data buffer 9-dimensional
    allocate(buff_send9_j(xs(1),xs(3),9,2*nHalo), &
             buff_send9_k(xs(1),xs(2),9,2*nHalo))
    ! receiving data buffer 9-dimensional
    allocate(buff_recv9_j(xs(1),xs(3),9,2*nHalo), &
             buff_recv9_k(xs(1),xs(2),9,2*nHalo))
    ! sending data buffer 1-dimensional
    allocate(buff_send1_j(xs(1),xs(3),2*nHalo), &
             buff_send1_k(xs(1),xs(2),2*nHalo))
    ! receiving data buffer 1-dimensional
    allocate(buff_recv1_j(xs(1),xs(3),2*nHalo), &
             buff_recv1_k(xs(1),xs(2),2*nHalo))
#else    
    ! no data buffer for CPU
    allocate(senjP(xs(1),xs(3),nHalo), &
             rcvjP(xs(1),xs(3),nHalo), &
             senjM(xs(1),xs(3),nHalo), &
             rcvjM(xs(1),xs(3),nHalo), &
             senkP(xs(1),xs(2),nHalo), &
             rcvkP(xs(1),xs(2),nHalo), &
             senkM(xs(1),xs(2),nHalo), &
             rcvkM(xs(1),xs(2),nHalo)  )
#endif    
    ! p_row: number of computational partition units in spanwise y-direction
    ! p_col: number of computational partition units in streamise z-direction
    dimens(1) = p_row
    dimens(2) = p_col
    neigh%inlet  = .false.
    neigh%outlet = .false.
    ! new communicator for cartesian topology
    call mpi_cart_coords(comm_cart,myid,2,coordr,ierr)
    coordn(1)=coordr(1)
    coords(1)=coordr(1)
    coordt(2)=coordr(2)
    coordb(2)=coordr(2)
    coordn(2)=coordr(2)+1
    if (coordn(2).gt.p_col-1) then
      coordn(2) = 0
      neigh%outlet = .true.
    endif
    coords(2)=coordr(2)-1
    if (coords(2).lt.0) then
      coords(2) = p_col-1
      neigh%inlet = .true.
    endif
    coordt(1)=coordr(1)+1
    if (coordt(1).gt.p_row-1) coordt(1) = 0
    coordb(1)=coordr(1)-1
    if (coordb(1).lt.0) coordb(1) = p_row-1
    neigh%myid = myid
    ! determine process ranks in the two MPI directions (y and z)
    call mpi_cart_rank(comm_cart, coordn, neigh%kp, ierr)
    call mpi_cart_rank(comm_cart, coords, neigh%km, ierr)
    call mpi_cart_rank(comm_cart, coordt, neigh%jp, ierr)
    call mpi_cart_rank(comm_cart, coordb, neigh%jm, ierr)
    ! different treatment of the spanwise direction for 2-D and 3-D simulations
    if (xs(2)>1) then
      if ((p_row>1) .and. (xs(2)>6)) then
        char_buff_j='mpi_copy_3d'
#define ACC_MPI_COPY_3D
      else if ((p_row==1) .and. (xs(2)>6)) then
        char_buff_j='local_copy_3d'
#define ACC_LOCAL_COPY_3D
      else if ((p_row==1) .and. (xs(2)<6)) then
        if (myid == 0) write(*,*) "ERROR: too few points jmax<2 x Halo, increase jmax!"
          call decomp_2d_finalize
          call mpi_finalize(ierr) 
          stop 
      endif
    else if (xs(2)==1) then
      char_buff_j='local_copy_2d'
#define ACC_LOCAL_COPY_2D
    endif
    ! size check
    if ((jmax<p_row) .or. (kmax<p_col)) then
      if ((myid == 0) .and. (jmax<p_row)) write(*,*) "ERROR: n_row>jmax, increase the number of points in y!"
      if ((myid == 0) .and. (kmax<p_col)) write(*,*) "ERROR: n_col>kmax, increase the number of points in z!"
      call decomp_2d_finalize
      call mpi_finalize(ierr) 
      stop 
    endif
! define GPU task IDs and initialize OpenACC runtime
#if defined(_OPENACC)
    block
      use openacc
      integer :: local_comm,mydev,ierr,ndev
      integer(acc_device_kind) :: dev_type
      call MPI_COMM_SPLIT_TYPE(comm_cart,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,local_comm,ierr)
      call MPI_COMM_RANK(local_comm,mydev,ierr)
#if defined (_NVIDIA)
      dev_type = acc_get_device_type()
      ndev  = acc_get_num_devices(dev_type)
      mydev = mod(mydev,ndev)
      call acc_set_device_num(mydev,dev_type)
      call acc_init(dev_type)
    end block
#endif
#if defined (_AMD)
      mydev = mydev
      call acc_init(dev_type)
      call acc_set_device_num(mydev,dev_type)
    end block
#endif
#endif
  !$acc enter data copyin(xs,char_buff_j) async
  !$acc enter data create(buff_send9_j,buff_send9_k,buff_recv9_j,buff_recv9_k) 
  !$acc enter data create(buff_send1_j,buff_send1_k,buff_recv1_j,buff_recv1_k)
  end subroutine
! initialization of the rescale MPI routine
  subroutine comm_init_rescale(z,xs)
  use decomp_2d
  use mod_param
  implicit none
  integer               :: i,k,ierr,rcy_col
  integer, dimension(2) :: coordr,coord_rcy,coord_inl
  integer, dimension(3) :: xs
  real(mytype), dimension(:) :: z
  ! sending data buffer 1-dimensional
  allocate(send_rcy(xs(1),xs(2),5,nHalo+1))
  ! receiving data buffer 1-dimensional
  allocate(recv_rcy(xs(1),xs(2),5,nHalo+1))
  ! checking whether the recycle length is appropriate
  if (pert_calc==1) then
    write(stdout,* ) 'ERROR: Please turn off perturbations when you use rescale-reintroduction'
    call decomp_2d_finalize
    call mpi_finalize(ierr) 
    stop
  endif
  if (spInlLen .gt. 0.0_mytype) then
    write(stdout,* ) 'ERROR: Please turn off sponge at inlet when you use rescale-reintroduction'
    call decomp_2d_finalize
    call mpi_finalize(ierr) 
    stop
  endif
  if ((z_recycle .lt. spInlLen) .or. (z_recycle .gt. len_z-spOutLen))then
    write(stdout,* ) 'ERROR: Rescycle location is inside the sponge'
    call decomp_2d_finalize
    call mpi_finalize(ierr) 
    stop
  endif
  ! Detect rank and k-index of the recycle postion
  neigh%recycle = .false.
  rcy_col = 0
  krcy  = 0
  if ((z(1) .le. z_recycle) .and. (z(xs(3)) .ge. z_recycle))then
    neigh%recycle = .true.
    call mpi_cart_coords(DECOMP_2D_COMM_CART_X,nrank,2,coordr,ierr)
    rcy_col = coordr(2)
    do k=1,xs(3)
      if (z(k) .gt. z_recycle) then
        krcy = k
        exit
      endif
    enddo
  endif
  call mpi_allreduce(MPI_IN_PLACE, rcy_col, 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr);   rcy_col = rcy_col/p_row
  call mpi_allreduce(MPI_IN_PLACE, krcy   , 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr);   krcy    = krcy   /p_row
  ! Set MPI communications
  if ((neigh%inlet) .or. (neigh%recycle)) then
    call mpi_cart_coords(DECOMP_2D_COMM_CART_X,nrank,2,coordr,ierr)
    coord_rcy(1)=coordr(1)
    coord_inl(1)=coordr(1)
    coord_rcy(2)=coordr(2)+rcy_col
    if (coord_rcy(2).gt.rcy_col) then
      coord_rcy(2) = 0
    endif
    coord_inl(2)=coordr(2)-rcy_col
    if (coord_inl(2).lt.0) then
      coord_inl(2) = rcy_col
    endif
    call mpi_cart_rank(DECOMP_2D_COMM_CART_X, coord_rcy, neigh%kp_rcy, ierr)
    call mpi_cart_rank(DECOMP_2D_COMM_CART_X, coord_inl, neigh%km_rcy, ierr)
  endif 
  !$acc enter data create(send_rcy,recv_rcy) 
end subroutine
! update halo cells in j-and k-direction only for GPU architecture   
#if defined(_OPENACC)
  ! 9-dimensional buffer
  subroutine haloUpdate9_jk(dir,xs,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9)
    use mod_param
    implicit none
    integer :: i,j,k,h,ierr,istat(mpi_status_size)
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: &
            tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
    logical, dimension(3) :: dir
    integer, dimension(3) :: xs
    ! sending or copying the jM & jP depending on *case*
    select case (char_buff_j)
      case("mpi_copy_3d")
        !$acc parallel loop collapse(3) default(present) async(2)
        do k=1,xs(3)
          do h=1,nHalo
            do i=1,xs(1)
              ! sendjM 
              buff_send9_j(i,k,1,h)       = tmp1(i,xs(2)+1-h,k) 
              buff_send9_j(i,k,2,h)       = tmp2(i,xs(2)+1-h,k) 
              buff_send9_j(i,k,3,h)       = tmp3(i,xs(2)+1-h,k)  
              buff_send9_j(i,k,4,h)       = tmp4(i,xs(2)+1-h,k) 
              buff_send9_j(i,k,5,h)       = tmp5(i,xs(2)+1-h,k)
              buff_send9_j(i,k,6,h)       = tmp6(i,xs(2)+1-h,k)  
              buff_send9_j(i,k,7,h)       = tmp7(i,xs(2)+1-h,k) 
              buff_send9_j(i,k,8,h)       = tmp8(i,xs(2)+1-h,k) 
              buff_send9_j(i,k,9,h)       = tmp9(i,xs(2)+1-h,k) 
              ! sendjP
              buff_send9_j(i,k,1,h+nHalo) = tmp1(i,h,k)         
              buff_send9_j(i,k,2,h+nHalo) = tmp2(i,h,k)              
              buff_send9_j(i,k,3,h+nHalo) = tmp3(i,h,k)                
              buff_send9_j(i,k,4,h+nHalo) = tmp4(i,h,k)                 
              buff_send9_j(i,k,5,h+nHalo) = tmp5(i,h,k)                   
              buff_send9_j(i,k,6,h+nHalo) = tmp6(i,h,k)                     
              buff_send9_j(i,k,7,h+nHalo) = tmp7(i,h,k)                    
              buff_send9_j(i,k,8,h+nHalo) = tmp8(i,h,k)                   
              buff_send9_j(i,k,9,h+nHalo) = tmp9(i,h,k)                                 
            enddo
          enddo
        enddo
      case("local_copy_3d")
        !$acc parallel loop collapse(3) default(present) async(2)
        do k=1,xs(3)
          do h=1,nHalo
            do i=1,xs(1)
              ! copyjM 
              tmp1(i,xs(2)+h,k) = tmp1(i,h,k)
              tmp2(i,xs(2)+h,k) = tmp2(i,h,k)
              tmp3(i,xs(2)+h,k) = tmp3(i,h,k)
              tmp4(i,xs(2)+h,k) = tmp4(i,h,k)
              tmp5(i,xs(2)+h,k) = tmp5(i,h,k)
              tmp6(i,xs(2)+h,k) = tmp6(i,h,k)
              tmp7(i,xs(2)+h,k) = tmp7(i,h,k)
              tmp8(i,xs(2)+h,k) = tmp8(i,h,k)
              tmp9(i,xs(2)+h,k) = tmp9(i,h,k)
              ! copyjP 
              tmp1(i,1-h,k) = tmp1(i,xs(2)+1-h,k) 
              tmp2(i,1-h,k) = tmp2(i,xs(2)+1-h,k)            
              tmp3(i,1-h,k) = tmp3(i,xs(2)+1-h,k)              
              tmp4(i,1-h,k) = tmp4(i,xs(2)+1-h,k)              
              tmp5(i,1-h,k) = tmp5(i,xs(2)+1-h,k)              
              tmp6(i,1-h,k) = tmp6(i,xs(2)+1-h,k)             
              tmp7(i,1-h,k) = tmp7(i,xs(2)+1-h,k)             
              tmp8(i,1-h,k) = tmp8(i,xs(2)+1-h,k)              
              tmp9(i,1-h,k) = tmp9(i,xs(2)+1-h,k)
           enddo
          enddo
        enddo
      case("local_copy_2d")
        !$acc parallel loop collapse(3) default(present) async(2)
        do k=1,xs(3)
          do h=1,nHalo
            do i=1,xs(1)
              ! copyjM 
              tmp1(i,1-h,k) = tmp1(i,1,k)
              tmp2(i,1-h,k) = tmp2(i,1,k)
              tmp3(i,1-h,k) = tmp3(i,1,k)
              tmp4(i,1-h,k) = tmp4(i,1,k)
              tmp5(i,1-h,k) = tmp5(i,1,k)
              tmp6(i,1-h,k) = tmp6(i,1,k)
              tmp7(i,1-h,k) = tmp7(i,1,k)
              tmp8(i,1-h,k) = tmp8(i,1,k)
              tmp9(i,1-h,k) = tmp9(i,1,k)
              ! copyjP
              tmp1(i,xs(2)+h,k) = tmp1(i,1,k)              
              tmp2(i,xs(2)+h,k) = tmp2(i,1,k)
              tmp3(i,xs(2)+h,k) = tmp3(i,1,k)             
              tmp4(i,xs(2)+h,k) = tmp4(i,1,k)            
              tmp5(i,xs(2)+h,k) = tmp5(i,1,k)            
              tmp6(i,xs(2)+h,k) = tmp6(i,1,k)          
              tmp7(i,xs(2)+h,k) = tmp7(i,1,k)          
              tmp8(i,xs(2)+h,k) = tmp8(i,1,k)            
              tmp9(i,xs(2)+h,k) = tmp9(i,1,k)
           enddo
          enddo
        enddo
    end select 
    ! sending the kM & kP
    !$acc parallel loop collapse(3) default(present) async(2)
    do j=1,xs(2)
      do h=1,nHalo
        do i=1,xs(1)
          ! sendkM 
          buff_send9_k(i,j,1,h)       = tmp1(i,j,xs(3)+1-h) 
          buff_send9_k(i,j,2,h)       = tmp2(i,j,xs(3)+1-h)
          buff_send9_k(i,j,3,h)       = tmp3(i,j,xs(3)+1-h)  
          buff_send9_k(i,j,4,h)       = tmp4(i,j,xs(3)+1-h) 
          buff_send9_k(i,j,5,h)       = tmp5(i,j,xs(3)+1-h)  
          buff_send9_k(i,j,6,h)       = tmp6(i,j,xs(3)+1-h) 
          buff_send9_k(i,j,7,h)       = tmp7(i,j,xs(3)+1-h)  
          buff_send9_k(i,j,8,h)       = tmp8(i,j,xs(3)+1-h)  
          buff_send9_k(i,j,9,h)       = tmp9(i,j,xs(3)+1-h)
          ! sendkP
          buff_send9_k(i,j,1,h+nHalo) = tmp1(i,j,h)         
          buff_send9_k(i,j,2,h+nHalo) = tmp2(i,j,h)                 
          buff_send9_k(i,j,3,h+nHalo) = tmp3(i,j,h)                 
          buff_send9_k(i,j,4,h+nHalo) = tmp4(i,j,h)                 
          buff_send9_k(i,j,5,h+nHalo) = tmp5(i,j,h)                  
          buff_send9_k(i,j,6,h+nHalo) = tmp6(i,j,h)                      
          buff_send9_k(i,j,7,h+nHalo) = tmp7(i,j,h)                
          buff_send9_k(i,j,8,h+nHalo) = tmp8(i,j,h)                     
          buff_send9_k(i,j,9,h+nHalo) = tmp9(i,j,h)         
        enddo
      enddo
    enddo
    !$acc wait(2)
! MPI call
#if defined(_GPU_DIRECT) && defined(ACC_MPI_COPY_3D)
  !$acc host_data use_device(buff_send9_k,buff_send9_j,buff_recv9_j,buff_recv9_k)
#elif defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_3D)
  !$acc host_data use_device(buff_send9_k,buff_recv9_k)
#elif defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_2D)
  !$acc host_data use_device(buff_send9_k,buff_recv9_k)
#elif !defined(_GPU_DIRECT) && defined(ACC_MPI_COPY_3D)  
  !$acc update host(buff_send9_k,buff_send9_j)
#elif !defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_3D)
  !$acc update host(buff_send9_k)  
#elif !defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_2D)
  !$acc update host(buff_send9_k)
#endif
    select case (char_buff_j)
      case("mpi_copy_3d")  
        ! MPI jM
        call mpi_sendrecv(buff_send9_j(:,:,:,1:nHalo),xs(1)*xs(3)*9*nHalo,real_type,neigh%jp,0, &
                          buff_recv9_j(:,:,:,1:nHalo),xs(1)*xs(3)*9*nHalo,real_type,neigh%jm,0, &
                          MPI_COMM_WORLD,istat,ierr) 
        ! MPI jP
        call mpi_sendrecv(buff_send9_j(:,:,:,nHalo+1:nHalo*2),xs(1)*xs(3)*9*nHalo,real_type,neigh%jm,0, &
                          buff_recv9_j(:,:,:,nHalo+1:nHalo*2),xs(1)*xs(3)*9*nHalo,real_type,neigh%jp,0, &
                          MPI_COMM_WORLD,istat,ierr) 
    end select
    ! MPI kM
    call mpi_sendrecv(buff_send9_k(:,:,:,1:nHalo),xs(1)*xs(2)*9*nHalo,real_type,neigh%kp,0, &
                      buff_recv9_k(:,:,:,1:nHalo),xs(1)*xs(2)*9*nHalo,real_type,neigh%km,0, &
                      MPI_COMM_WORLD,istat,ierr) 
    ! MPI kP
    call mpi_sendrecv(buff_send9_k(:,:,:,nHalo+1:2*nHalo),xs(1)*xs(2)*9*nHalo,real_type,neigh%km,0, &
                      buff_recv9_k(:,:,:,nHalo+1:2*nHalo),xs(1)*xs(2)*9*nHalo,real_type,neigh%kp,0, &
                      MPI_COMM_WORLD,istat,ierr) 
! MPI call
#if defined(_GPU_DIRECT)
      !$acc end host_data
#elif !defined(_GPU_DIRECT) && defined(ACC_MPI_COPY_3D)
      !$acc update device(buff_recv9_k,buff_recv9_j)
#elif !defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_3D)
      !$acc update device(buff_recv9_k)
#elif !defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_2D)
      !$acc update device(buff_recv9_k)
#endif
!$acc wait(2)
    ! receving the jM & jP depending on *case*
    select case (char_buff_j)
      case("mpi_copy_3d")
        !$acc parallel loop collapse(3) default(present) async(2)
        do k=1,xs(3)
          do h=1,nHalo
            do i=1,xs(1)
              ! rcvjM
              tmp1(i,1-h,k) = buff_recv9_j(i,k,1,h) 
              tmp2(i,1-h,k) = buff_recv9_j(i,k,2,h)
              tmp3(i,1-h,k) = buff_recv9_j(i,k,3,h) 
              tmp4(i,1-h,k) = buff_recv9_j(i,k,4,h) 
              tmp5(i,1-h,k) = buff_recv9_j(i,k,5,h)
              tmp6(i,1-h,k) = buff_recv9_j(i,k,6,h)  
              tmp7(i,1-h,k) = buff_recv9_j(i,k,7,h) 
              tmp8(i,1-h,k) = buff_recv9_j(i,k,8,h)
              tmp9(i,1-h,k) = buff_recv9_j(i,k,9,h) 
              ! rcvjP
              tmp1(i,xs(2)+h,k) = buff_recv9_j(i,k,1,h+nHalo) 
              tmp2(i,xs(2)+h,k) = buff_recv9_j(i,k,2,h+nHalo) 
              tmp3(i,xs(2)+h,k) = buff_recv9_j(i,k,3,h+nHalo) 
              tmp4(i,xs(2)+h,k) = buff_recv9_j(i,k,4,h+nHalo) 
              tmp5(i,xs(2)+h,k) = buff_recv9_j(i,k,5,h+nHalo) 
              tmp6(i,xs(2)+h,k) = buff_recv9_j(i,k,6,h+nHalo)             
              tmp7(i,xs(2)+h,k) = buff_recv9_j(i,k,7,h+nHalo)              
              tmp8(i,xs(2)+h,k) = buff_recv9_j(i,k,8,h+nHalo)              
              tmp9(i,xs(2)+h,k) = buff_recv9_j(i,k,9,h+nHalo) 
            enddo
          enddo
        enddo
    end select  
    ! receving the kM & kP 
    !$acc parallel loop collapse(3) default(present) async(2)
    do j=1,xs(2)
      do h=1,nHalo
        do i=1,xs(1)
          ! rcvkM
          tmp1(i,j,1-h) = buff_recv9_k(i,j,1,h) 
          tmp2(i,j,1-h) = buff_recv9_k(i,j,2,h)
          tmp3(i,j,1-h) = buff_recv9_k(i,j,3,h) 
          tmp4(i,j,1-h) = buff_recv9_k(i,j,4,h)
          tmp5(i,j,1-h) = buff_recv9_k(i,j,5,h)
          tmp6(i,j,1-h) = buff_recv9_k(i,j,6,h) 
          tmp7(i,j,1-h) = buff_recv9_k(i,j,7,h) 
          tmp8(i,j,1-h) = buff_recv9_k(i,j,8,h)
          tmp9(i,j,1-h) = buff_recv9_k(i,j,9,h)  
          ! rcvkP
          tmp1(i,j,xs(3)+h) = buff_recv9_k(i,j,1,h+nHalo) 
          tmp2(i,j,xs(3)+h) = buff_recv9_k(i,j,2,h+nHalo)      
          tmp3(i,j,xs(3)+h) = buff_recv9_k(i,j,3,h+nHalo)        
          tmp4(i,j,xs(3)+h) = buff_recv9_k(i,j,4,h+nHalo)          
          tmp5(i,j,xs(3)+h) = buff_recv9_k(i,j,5,h+nHalo)         
          tmp6(i,j,xs(3)+h) = buff_recv9_k(i,j,6,h+nHalo)         
          tmp7(i,j,xs(3)+h) = buff_recv9_k(i,j,7,h+nHalo)         
          tmp8(i,j,xs(3)+h) = buff_recv9_k(i,j,8,h+nHalo)          
          tmp9(i,j,xs(3)+h) = buff_recv9_k(i,j,9,h+nHalo) 
        enddo
      enddo
    enddo
  end subroutine
  ! 1-dimensional buffer
  subroutine haloUpdate1_jk(dir,xs,tmp1)
    use mod_param
    implicit none
    integer :: i,j,k,h,ierr,istat(mpi_status_size)
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: tmp1
    logical, dimension(3) :: dir
    integer, dimension(3) :: xs
    ! sending or copying the jM & jP depending on *case*
    select case (char_buff_j)
      case("mpi_copy_3d")
        !$acc parallel loop collapse(3) default(present) async(2)
        do k=1,xs(3)
          do h=1,nHalo
            do i=1,xs(1)
              ! sendjM
              buff_send1_j(i,k,h)       = tmp1(i,xs(2)+1-h,k) 
              ! sendjP 
              buff_send1_j(i,k,h+nHalo) = tmp1(i,h,k)        
            enddo
          enddo
        enddo
      case("local_copy_3d")
        !$acc parallel loop collapse(3) default(present) async(2)
        do k=1,xs(3)
          do h=1,nHalo
            do i=1,xs(1)
              ! copyjM
              tmp1(i,xs(2)+h,k) = tmp1(i,h,k)
              ! copyjP
              tmp1(i,1-h,k) = tmp1(i,xs(2)+1-h,k)
            enddo
          enddo
        enddo
      case("local_copy_2d")
        !$acc parallel loop collapse(3) default(present) async(2)
        do k=1,xs(3)
          do h=1,nHalo
            do i=1,xs(1)
              ! copyjM
              tmp1(i,1-h,k) = tmp1(i,1,k)
              ! copyjP
              tmp1(i,xs(2)+h,k) = tmp1(i,1,k)
            enddo
          enddo
        enddo
    end select
    ! sending the kM & kP
    !$acc parallel loop collapse(3) default(present) async(2)
    do j=1,xs(2)
      do h=1,nHalo
        do i=1,xs(1)
          buff_send1_k(i,j,h)       = tmp1(i,j,xs(3)+1-h) 
          buff_send1_k(i,j,h+nHalo) = tmp1(i,j,h)        
        enddo
      enddo
    enddo
    !$acc wait(2)
! MPI call
#if defined(_GPU_DIRECT) && defined(ACC_MPI_COPY_3D)
  !$acc host_data use_device(buff_send1_k,buff_send1_j,buff_recv1_j,buff_recv1_k)
#elif defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_3D)
  !$acc host_data use_device(buff_send1_k,buff_recv1_k)
#elif defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_2D)
  !$acc host_data use_device(buff_send1_k,buff_recv1_k)
#elif !defined(_GPU_DIRECT) && defined(ACC_MPI_COPY_3D)  
  !$acc update host(buff_send1_k,buff_send1_j)
#elif !defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_3D)
  !$acc update host(buff_send1_k)  
#elif !defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_2D)
  !$acc update host(buff_send1_k)
#endif
    select case (char_buff_j)
      case("mpi_copy_3d")
        ! MPI jM
        call mpi_sendrecv(buff_send1_j(:,:,1:nHalo),xs(1)*xs(3)*nHalo,real_type,neigh%jp,0, &
                          buff_recv1_j(:,:,1:nHalo),xs(1)*xs(3)*nHalo,real_type,neigh%jm,0, &
                          MPI_COMM_WORLD,istat,ierr)
        ! MPI jP
        call mpi_sendrecv(buff_send1_j(:,:,nHalo+1:nHalo*2),xs(1)*xs(3)*nHalo,real_type,neigh%jm,0, &
                          buff_recv1_j(:,:,nHalo+1:nHalo*2),xs(1)*xs(3)*nHalo,real_type,neigh%jp,0, &
                          MPI_COMM_WORLD,istat,ierr) 
    end select 
    ! MPI kM 
    call mpi_sendrecv(buff_send1_k(:,:,1:nHalo),xs(1)*xs(2)*nHalo,real_type,neigh%kp,0, &
                      buff_recv1_k(:,:,1:nHalo),xs(1)*xs(2)*nHalo,real_type,neigh%km,0, &
                      MPI_COMM_WORLD,istat,ierr) 
    ! MPI kP
    call mpi_sendrecv(buff_send1_k(:,:,nHalo+1:2*nHalo),xs(1)*xs(2)*nHalo,real_type,neigh%km,0, &
                      buff_recv1_k(:,:,nHalo+1:2*nHalo),xs(1)*xs(2)*nHalo,real_type,neigh%kp,0, &
                      MPI_COMM_WORLD,istat,ierr) 
#if defined(_GPU_DIRECT)
      !$acc end host_data
#elif !defined(_GPU_DIRECT) && defined(ACC_MPI_COPY_3D)
      !$acc update device(buff_recv1_k,buff_recv1_j)
#elif !defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_3D)
      !$acc update device(buff_recv1_k)
#elif !defined(_GPU_DIRECT) && defined(ACC_LOCAL_COPY_2D)
      !$acc update device(buff_recv1_k)
#endif
!$acc wait(2)
    ! receving the jM & jP depending on *case*
    select case (char_buff_j)
      case("mpi_copy_3d")
        !$acc parallel loop collapse(3) default(present) async(1)
        do k=1,xs(3)
          do h=1,nHalo
            do i=1,xs(1)
              ! rcvjM
              tmp1(i,1-h,k) = buff_recv1_j(i,k,h) 
              ! rcvjP
              tmp1(i,xs(2)+h,k) = buff_recv1_j(i,k,h+nHalo) 
            enddo
          enddo
        enddo
    end select
    ! receving the kM & kP 
    !$acc parallel loop collapse(3) default(present) async(1)
    do j=1,xs(2)
      do h=1,nHalo
        do i=1,xs(1)
          ! rcvkM
          tmp1(i,j,1-h) = buff_recv1_k(i,j,h)
          ! rcvkP
          tmp1(i,j,xs(3)+h) = buff_recv1_k(i,j,h+nHalo)
        enddo
      enddo
    enddo
  end subroutine
! update halo cells in i-direction only for GPU architecture 
  subroutine haloUpdateMult_i(dir,xs,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9)
    use mod_param
    implicit none
    integer :: h,i,j,k
    logical, dimension(3) :: dir
    integer, dimension(3) :: xs
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:), optional :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
    ! copying the iM & iP depending on how many tmp present
    if (present(tmp1)) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,xs(3)
        do j=1,xs(2)   
          do h=1,nHalo
            tmp1(xs(1)+h,j,k) = tmp1(h,j,k)
            tmp1(1-h,j,k) = tmp1(xs(1)+1-h,j,k)
          enddo 
        enddo
      enddo
    endif
    if (present(tmp2)) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,xs(3)
        do j=1,xs(2)   
          do h=1,nHalo
            tmp2(xs(1)+h,j,k) = tmp2(h,j,k)
            tmp2(1-h,j,k) = tmp2(xs(1)+1-h,j,k)
          enddo 
        enddo
      enddo
    endif
    if (present(tmp3)) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,xs(3)
        do j=1,xs(2)   
          do h=1,nHalo
            tmp3(xs(1)+h,j,k) = tmp3(h,j,k)
            tmp3(1-h,j,k) = tmp3(xs(1)+1-h,j,k)
          enddo 
        enddo
      enddo
    endif
    if (present(tmp4)) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,xs(3)
        do j=1,xs(2)   
          do h=1,nHalo
            tmp4(xs(1)+h,j,k) = tmp4(h,j,k)
            tmp4(1-h,j,k) = tmp4(xs(1)+1-h,j,k)
          enddo 
        enddo
      enddo
    endif
    if (present(tmp5)) then
      !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,xs(3)
      do j=1,xs(2)   
        do h=1,nHalo
          tmp5(xs(1)+h,j,k) = tmp5(h,j,k)
          tmp5(1-h,j,k) = tmp5(xs(1)+1-h,j,k)
        enddo 
      enddo
    enddo
    endif
    if (present(tmp6)) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,xs(3)
        do j=1,xs(2)   
          do h=1,nHalo
            tmp6(xs(1)+h,j,k) = tmp6(h,j,k)
            tmp6(1-h,j,k) = tmp6(xs(1)+1-h,j,k)
          enddo 
        enddo
      enddo
    endif
    if (present(tmp7)) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,xs(3)
        do j=1,xs(2)   
          do h=1,nHalo
            tmp7(xs(1)+h,j,k) = tmp7(h,j,k)
            tmp7(1-h,j,k) = tmp7(xs(1)+1-h,j,k)
          enddo 
        enddo
      enddo
    endif
    if (present(tmp8)) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,xs(3)
        do j=1,xs(2)   
          do h=1,nHalo
            tmp8(xs(1)+h,j,k) = tmp8(h,j,k)
            tmp8(1-h,j,k) = tmp8(xs(1)+1-h,j,k)
          enddo 
        enddo
      enddo
    endif
    if (present(tmp9)) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,xs(3)
        do j=1,xs(2)   
          do h=1,nHalo
            tmp9(xs(1)+h,j,k) = tmp9(h,j,k)
            tmp9(1-h,j,k) = tmp9(xs(1)+1-h,j,k)
          enddo 
        enddo
      enddo
    endif 
  end subroutine
! update inlet location with recycling position
  subroutine haloUpdate_rescale(xs,tmp1,tmp2,tmp3,tmp4,tmp5)
    use mod_param
    implicit none
    integer :: ierr,istat(mpi_status_size),i,j,h
    integer, dimension(3) :: xs
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: tmp1,tmp2,tmp3,tmp4,tmp5
    ! sending variables from recycle location to inlet
    !$acc parallel loop collapse(3) default(present)
    do j=1,xs(2)
      do h=1,nHalo+1
        do i=1,xs(1)
          send_rcy(i,j,1,h) = tmp1(i,j,krcy+1-h)
          send_rcy(i,j,2,h) = tmp2(i,j,krcy+1-h)
          send_rcy(i,j,3,h) = tmp3(i,j,krcy+1-h)
          send_rcy(i,j,4,h) = tmp4(i,j,krcy+1-h)
          send_rcy(i,j,5,h) = tmp5(i,j,krcy+1-h)
        enddo
      enddo
    enddo
! MPI call
#if defined(_GPU_DIRECT)
  !$acc host_data use_device(send_rcy,recv_rcy)  
#elif !defined(_GPU_DIRECT)
  !$acc update host(send_rcy)
#endif
    call mpi_sendrecv(send_rcy,xs(1)*xs(2)*5*(nHalo+1),real_type,neigh%km_rcy,0, &
                      recv_rcy,xs(1)*xs(2)*5*(nHalo+1),real_type,neigh%kp_rcy,0, &
                      MPI_COMM_WORLD,istat,ierr)
#if defined(_GPU_DIRECT)
  !$acc end host_data  
#elif !defined(_GPU_DIRECT)
  !$acc update device(recv_rcy) 
#endif
    ! updating halo and inlet
    if((neigh%inlet)) then
      !$acc parallel loop collapse(3) default(present)
      do j=1,xs(2)
        do h=1,nHalo+1
          do i=1,xs(1)
            tmp1(i,j,2-h) = recv_rcy(i,j,1,h)
            tmp2(i,j,2-h) = recv_rcy(i,j,2,h)
            tmp3(i,j,2-h) = recv_rcy(i,j,3,h)
            tmp4(i,j,2-h) = recv_rcy(i,j,4,h)
            tmp5(i,j,2-h) = recv_rcy(i,j,5,h)
          enddo 
        enddo
      enddo
    endif
  end subroutine
! update halo cells only for CPU architecture   
#else
  subroutine haloUpdateMult_CPU(dir,xs,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9)
    use mod_param
    implicit none
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:), optional :: &
            tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
    logical, dimension(3) :: dir
    integer, dimension(3) :: xs
    ! call subroutine haloUpdate
    if(present(tmp1)) call haloUpdate(dir,xs,tmp1)
    if(present(tmp2)) call haloUpdate(dir,xs,tmp2)
    if(present(tmp3)) call haloUpdate(dir,xs,tmp3)
    if(present(tmp4)) call haloUpdate(dir,xs,tmp4)
    if(present(tmp5)) call haloUpdate(dir,xs,tmp5)
    if(present(tmp6)) call haloUpdate(dir,xs,tmp6)
    if(present(tmp7)) call haloUpdate(dir,xs,tmp7)
    if(present(tmp8)) call haloUpdate(dir,xs,tmp8)
    if(present(tmp9)) call haloUpdate(dir,xs,tmp9)
  end subroutine
  ! subroutine valid for all three directions
  subroutine haloUpdate(dir,xs,array)
    use mod_param
    implicit none
    integer :: ierr,istat(mpi_status_size), i,j,k,h
    logical, dimension(3) :: dir
    integer, dimension(3) :: xs
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: array
    ! copy iM & iP
    if (dir(1)) then
      do h=1,nHalo
        array(xs(1)+h,1:xs(2),1:xs(3)) = array(           h,1:xs(2),1:xs(3))
        array(1-h,       1:xs(2),1:xs(3)) = array(xs(1)+1-h,1:xs(2),1:xs(3))
      enddo
    endif
    ! sending the jM & jP, and receiving them
    if (dir(2)) then 
      if (xs(2)>1) then
        do h=1,nHalo
          senjM(1:xs(1),1:xs(3),h) = array(1:xs(1),xs(2)+1-h,1:xs(3))
          senjP(1:xs(1),1:xs(3),h) = array(1:xs(1),        h,1:xs(3))
        enddo
        call mpi_sendrecv(senjM,xs(1)*xs(3)*nHalo,real_type,neigh%jp,0, &
                          rcvjM,xs(1)*xs(3)*nHalo,real_type,neigh%jm,0, &
                          MPI_COMM_WORLD,istat,ierr)
        call mpi_sendrecv(senjP,xs(1)*xs(3)*nHalo,real_type,neigh%jm,0, &
                          rcvjP,xs(1)*xs(3)*nHalo,real_type,neigh%jp,0, &
                          MPI_COMM_WORLD,istat,ierr)
        do h=1,nHalo
          array(1:xs(1),1-h,    1:xs(3)) = rcvjM(1:xs(1),1:xs(3),h)
          array(1:xs(1),xs(2)+h,1:xs(3)) = rcvjP(1:xs(1),1:xs(3),h)
        enddo
      ! for 2D simulations just copy the data into the halo cells, then derivatives will be zero
      else 
        do h=1,nHalo 
          array(1:xs(1),1-h,    1:xs(3)) = array(1:xs(1),1,1:xs(3))
          array(1:xs(1),xs(2)+h,1:xs(3)) = array(1:xs(1),1,1:xs(3))
        enddo
      endif
    endif
    ! sending the kM & kP, and receiving them
    do h=1,nHalo 
      senkM(1:xs(1),1:xs(2),h) = array(1:xs(1),1:xs(2),xs(3)+1-h)
      senkP(1:xs(1),1:xs(2),h) = array(1:xs(1),1:xs(2),        h)
    enddo
    call mpi_sendrecv(senkM,xs(1)*xs(2)*nHalo,real_type,neigh%kp,0, &
                      rcvkM,xs(1)*xs(2)*nHalo,real_type,neigh%km,0, &
                      MPI_COMM_WORLD,istat,ierr)
    call mpi_sendrecv(senkP,xs(1)*xs(2)*nHalo,real_type,neigh%km,0, &
                      rcvkP,xs(1)*xs(2)*nHalo,real_type,neigh%kp,0, &
                      MPI_COMM_WORLD,istat,ierr)
    do h=1,nHalo
      array(1:xs(1),1:xs(2),1-h)     = rcvkM(1:xs(1),1:xs(2),h)
      array(1:xs(1),1:xs(2),xs(3)+h) = rcvkP(1:xs(1),1:xs(2),h)
    enddo
  end subroutine
! update inlet location with recycling position
  subroutine haloUpdate_rescale(xs,tmp1,tmp2,tmp3,tmp4,tmp5)
    use mod_param
    implicit none
    integer :: ierr,istat(mpi_status_size),h
    integer, dimension(3) :: xs
    real(mytype), dimension(1-nHalo:, 1-nHalo:, 1-nHalo:) :: tmp1,tmp2,tmp3,tmp4,tmp5
    ! sending variables from recycle location to inlet
    do h=1,nHalo+1
       send_rcy(1:xs(1),1:xs(2),1,h) = tmp1(1:xs(1),1:xs(2),krcy+1-h)
       send_rcy(1:xs(1),1:xs(2),2,h) = tmp2(1:xs(1),1:xs(2),krcy+1-h)
       send_rcy(1:xs(1),1:xs(2),3,h) = tmp3(1:xs(1),1:xs(2),krcy+1-h)
       send_rcy(1:xs(1),1:xs(2),4,h) = tmp4(1:xs(1),1:xs(2),krcy+1-h)
       send_rcy(1:xs(1),1:xs(2),5,h) = tmp5(1:xs(1),1:xs(2),krcy+1-h)
    enddo
    call mpi_sendrecv(send_rcy,xs(1)*xs(2)*5*(nHalo+1),real_type,neigh%km_rcy,0, &
                      recv_rcy,xs(1)*xs(2)*5*(nHalo+1),real_type,neigh%kp_rcy,0, &
                      MPI_COMM_WORLD,istat,ierr)
    ! updating halo and inlet
    if((neigh%inlet)) then
       do h=1,nHalo+1
          tmp1(1:xs(1),1:xs(2),2-h) = recv_rcy(1:xs(1),1:xs(2),1,h)
          tmp2(1:xs(1),1:xs(2),2-h) = recv_rcy(1:xs(1),1:xs(2),2,h)
          tmp3(1:xs(1),1:xs(2),2-h) = recv_rcy(1:xs(1),1:xs(2),3,h)
          tmp4(1:xs(1),1:xs(2),2-h) = recv_rcy(1:xs(1),1:xs(2),4,h)
          tmp5(1:xs(1),1:xs(2),2-h) = recv_rcy(1:xs(1),1:xs(2),5,h)
       enddo
    endif
  end subroutine
#endif
end module
