!----------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                                  !
!      Peng-Robinson routines                                                                                                      !
!                                                                                                                                  !
!      DESCRIPTION:                                                                                                                !
!      ------------                                                                                                                !
!      This module is for the Peng-Robinson EoS                                                                                    !
!                                                                                                                                  !
!                                                                                                                                  !
!      CHANGELOG:                                                                                                                  !
!      ----------                                                                                                                  !
!         xx.xx.2022: module created (Pietro)                                                                                      !
!                                                                                                                                  !
!----------------------------------------------------------------------------------------------------------------------------------!
  




  type(EOSModel_lookup) function init_EOSModel_lookup()
    init_EOSModel_lookup%name="look-up-table"
    call init_EOSModel_lookup%initEOS()
  end function init_EOSModel_lookup


  subroutine initialize_lookup(this)
    use mod_param
    use mod_table
    implicit none
    class(EOSModel_lookup) :: this

    call readInitTable_RE(RE_tab)

  end subroutine



  subroutine calcState_re_lookup(this,rho,ien,pre,tem,mu,ka,i1,i2,j1,j2,k1,k2)
    use decomp_2d
    use mpi
    use mod_param
    use mod_table
    implicit none
    class(EOSModel_lookup) :: this
    integer :: i,j,k,ierr
    integer, intent(IN) :: i1,i2,j1,j2,k1,k2
    real(mytype), intent(IN),  dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,ien
    real(mytype), intent(OUT), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: pre,tem,mu,ka
    real(mytype) :: err,errMax, errMin


!    !$acc parallel loop collapse(3) default(present) private(invrho)
    do k=k1,k2
      do j=j1,j2
        do i=i1,i2
          call tableLookUp(err, i,j,k, RE_tab, rho(i,j,k), ien(i,j,k), &
                            val1 = pre(i,j,k), val2 = tem(i,j,k), & 
                            val3 = mu(i,j,k),  val4 = ka(i,j,k)) 

          errMax = max(errMax,err)
          errMin = min(errMin,err)
        enddo
      enddo
    enddo

    call mpi_allreduce(errMax, err,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr); errMax = err
    call mpi_allreduce(errMin, err,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr); errMin = err

    if (errMax > 0) then 
      write(*,*) "ERROR CODE:  max/min ", errMax, errMin
      call decomp_2d_abort(-1, "Error in table interpolation!")
      stop
    endif

  end subroutine


  subroutine calcState_rP_lookup(this,rho,pre,ien,tem,mu,ka, i1,i2,j1,j2,k1,k2)
    use decomp_2d, only: mytype
    use mod_param
    use mod_table
    implicit none
    class(EOSModel_lookup) :: this
    integer :: i,j,k
    integer, intent(IN) :: i1,i2,j1,j2,k1,k2
    real(mytype), intent(IN),  dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,pre
    real(mytype), intent(OUT), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: ien,tem,mu,ka
    real(mytype) :: err

!    !$acc parallel loop collapse(3) default(present)
    do k=k1,k2
      do j=j1,j2
        do i=i1,i2
          call tableLookUp_Inverse(err, RE_tab, rho(i,j,k), RE_tab%val1, pre(i,j,k), ien(i,j,k))
        enddo 
      enddo 
    enddo 
  end subroutine


  subroutine calcState_rT_lookup(this,rho,tem,ien,pre,mu,ka, i1,i2,j1,j2,k1,k2)
    use decomp_2d, only: mytype
    use mod_param
    use mod_table
    implicit none
    class(EOSModel_lookup) :: this
    integer :: i,j,k
    integer, intent(IN) :: i1,i2,j1,j2,k1,k2
    real(mytype), intent(IN),  dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,tem
    real(mytype), intent(OUT), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: ien,pre,mu,ka
    real(mytype) :: err

!    !$acc parallel loop collapse(3) default(present)
    do k=k1,k2
      do j=j1,j2
        do i=i1,i2
          call tableLookUp_Inverse(err, RE_tab, rho(i,j,k), RE_tab%val2, tem(i,j,k), ien(i,j,k))
          call tableLookUp(err, i,j,k, RE_tab, rho(i,j,k), ien(i,j,k), &
                            val1=pre(i,j,k), val3=mu(i,j,k), val4=ka(i,j,k))
        enddo 
      enddo 
    enddo 
  end subroutine

  subroutine calcState_PT_lookup(this,pre,tem,rho,ien,mu,ka, i1,i2,j1,j2,k1,k2)
    use decomp_2d, only: mytype
    use mod_param
    use mod_table
    implicit none
    class(EOSModel_lookup) :: this
    integer :: i,j,k
    integer, intent(IN) :: i1,i2,j1,j2,k1,k2
    real(mytype), intent(IN),  dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: pre,tem
    real(mytype), intent(OUT), dimension(1-nHalo:,1-nHalo:,1-nHalo:) :: rho,ien,mu,ka
    real(mytype) :: err

    do k=k1,k2
      do j=j1,j2
        do i=i1,i2
          call tableLookUp(err, i,j,k, PT_tab, pre(i,j,k),   tem(i,j,k), val1=rho(i,j,k))
          call tableLookUp_Inverse(err, RE_tab, rho(i,j,k), RE_tab%val1, pre(i,j,k), ien(i,j,k))
          call tableLookUp(err, i,j,k, RE_tab, rho(i,j,k), ien(i,j,k), val3=mu(i,j,k), val4=ka(i,j,k))
        enddo 
      enddo 
    enddo 

  end subroutine


  subroutine calcSOS_re_lookup(this,rho,ien,sos)
    
    !$acc routine seq

    use decomp_2d, only: mytype
    use mod_param
    use mod_table

    implicit none
    class(EOSModel_lookup) :: this
    real(mytype), intent(IN)  :: rho,ien
    real(mytype), intent(OUT) :: sos
    real(mytype) :: err
    integer :: i,j,k
    i = 1
    j = 1
    k = 1

    call tableLookUp(err, i,j,k, RE_tab, rho, ien, val5=sos)

  end subroutine


  subroutine calcFac_re_lookup(this,rho,ien,fac)
    
    !$acc routine seq

    use decomp_2d, only: mytype
    use mod_param
    use mod_table

    implicit none
    class(EOSModel_lookup) :: this
    real(mytype), intent(IN)  :: rho,ien
    real(mytype), intent(OUT) :: fac
    real(mytype) :: err
    integer :: i,j,k
    i = 1
    j = 1
    k = 1
    

    call tableLookUp(err, i,j,k, RE_tab, rho, ien, val6=fac)

  end subroutine




