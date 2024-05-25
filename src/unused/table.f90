
module mod_table

  use decomp_2d
  use mod_param
  
  implicit none

  integer, parameter :: mytypeTable = 8

  type TABLE2D
    character(len=30) :: tableName
    integer :: diagonalTable
    integer :: imaxTab, jmaxTab, jmaxTabDiagonal
    real(mytypeTable), allocatable, dimension(:)   :: dim1, dim2
    real(mytypeTable), allocatable, dimension(:,:) :: val1,val2,val3,val4,val5,val6,val7,val8,val9,val10
  end type TABLE2D

  type (TABLE2D) :: RE_tab
  type (TABLE2D) :: PT_tab

contains




! READ AND INIT RHO-INT ENERGY TABLE
!
  subroutine readInitTable_RE(tab)
    use mod_param

    implicit none
    integer :: fs,i,j
    type (TABLE2D) :: tab

    tab%tableName = "RE table"

    call readTabData1D('tables/RE_tab_1d_rho.bin', tab%dim1, tab%imaxTab)
    call readTabData1D('tables/RE_tab_1d_ien.bin', tab%dim2, tab%jmaxTab)
    call readTabData2D('tables/RE_tab_2d_pre.bin', tab%val1, tab)
    call readTabData2D('tables/RE_tab_2d_tem.bin', tab%val2, tab)
    call readTabData2D('tables/RE_tab_2d_vis.bin', tab%val3, tab)
    call readTabData2D('tables/RE_tab_2d_con.bin', tab%val4, tab)
    call readTabData2D('tables/RE_tab_2d_sos.bin', tab%val5, tab)
    call readTabData2D('tables/RE_tab_2d_fac.bin', tab%val6, tab)
    call readTabData2D('tables/RE_tab_2d_der.bin', tab%val7, tab)
    
    ! call readTabData2D('tables/RE_tab_2d_dmudT.bin', tab%val7,  tab)
    ! call readTabData2D('tables/RE_tab_2d_dmudr.bin', tab%val8,  tab)
    ! call readTabData2D('tables/RE_tab_2d_dkadT.bin', tab%val9,  tab)
    ! call readTabData2D('tables/RE_tab_2d_dkadr.bin', tab%val10, tab)

  end subroutine


! READ AND INIT P-T TABLE
!
  subroutine readInitTable_PT(tab)
    use mod_param

    implicit none
    integer :: fs,i,j
    type (TABLE2D) :: tab

    tab%tableName = "PT table"

    call readTabData1D('tables/PT_tab_1d_pre.bin', tab%dim1, tab%imaxTab)
    call readTabData1D('tables/PT_tab_1d_tem.bin', tab%dim2, tab%jmaxTab)
    call readTabData2D('tables/PT_tab_2d_rho.bin', tab%val1, tab)

  end subroutine




  subroutine testTable(tab, val1_in, val2_in)

    use mod_param

    real(mytype)   :: rho, ien, error
    real(mytype)   :: val1_in, val2_in, val1, val2, val3, val4
    type (TABLE2D) :: tab

    call tableLookUp(error, 1,1,1,tab, val1_in, val2_in, val1=val1, val2=val2, val3=val3, val4=val4)

    if (nrank == 0) write(*,*) 'lookUpAllFromTable'
    if (nrank == 0) write(*,*) 'val1_in = ', val1_in, 'val2_in = ', val2_in
    if (nrank == 0) write(*,*) 'val1 = ', val1, 'val2 = ', val2
    if (nrank == 0) write(*,*) 'val3 = ', val3, 'val4 = ', val4

  end subroutine




  subroutine readTabData1D(name,tabData,ndata)
    implicit none
    logical :: file_exists
    integer :: n,ndata,fs,AllocateStatus,ierr
    character(len=*), intent(in) :: name
    real(mytypeTable), allocatable, dimension(:) :: tabData

    inquire(file = name, exist = file_exists)
    inquire(file = name, size = fs)
    ndata = fs/mytypeTable
    if (file_exists) then 
      allocate(tabData(ndata), STAT = AllocateStatus); if (AllocateStatus /= 0) STOP "Not enough memory for ien"
      open(11,file=name, access ='direct', recl=ndata*mytypeTable)
      read(11,rec=1) (tabData(n), n=1,ndata)
      close(11)
      if (nrank == 0) write(*,*) 'FINISHED READING: ', name, ', ndata = ', ndata, maxval(tabData), minval(tabData)
    else
      if (nrank == 0) write(*,*) 'FILE NOT FOUND: ', name
      call mpi_finalize(ierr)
      stop
    endif

  end subroutine


  subroutine readTabData2D(name,tabData,tab)
    implicit none
    logical :: file_exists
    integer :: i,j,AllocateStatus,ierr, fs,ndataPoints
    character(len=*), intent(in) :: name
    real(mytypeTable), allocatable, dimension(:,:) :: tabData
    type (TABLE2D) :: tab

    inquire(file = name, exist = file_exists)
    inquire(file = name, size = fs)
    ndataPoints = fs/mytypeTable

!   This is to check if the table has been stored in low storage diagonal format. If yes, then change the reading
!
    if (ndataPoints /= tab%imaxTab*tab%jmaxTab) then 
      tab%diagonalTable = 1
      tab%jmaxTabDiagonal = ndataPoints/tab%imaxTab
    else
      tab%diagonalTable = 0
      tab%jmaxTabDiagonal = tab%jmaxTab
    endif

    if (file_exists) then 
      allocate(tabData(tab%imaxTab,tab%jmaxTabDiagonal), STAT = AllocateStatus); if (AllocateStatus /= 0) STOP "Not enough memory"
      open(11,file=name, access ='direct', recl=tab%jmaxTabDiagonal*mytypeTable)
      do i=1,tab%imaxTab
        read(11,rec=i) (tabData(i,j), j=1,tab%jmaxTabDiagonal)
      enddo
      close(11)
      if (nrank == 0) write(*,*) 'FINISHED READING: ', name, ', size: ', tab%imaxTab, tab%jmaxTabDiagonal
    else
      if (nrank == 0) write(*,*) 'FILE NOT FOUND: ', name
      call mpi_finalize(ierr)
      stop
    endif

  end subroutine










  subroutine checkOutsideTable(iindex, jindex, tab, val1, val2, error, ii,jj,kk)

    implicit none
    integer :: iindex, jindex, jMin, jMax
    real(mytype) :: val1, val2
    real(mytype) :: error
    type (TABLE2D) :: tab
    integer :: ii,jj,kk

    error = 0

    if (iindex >= tab%imaxTab-1) then
      ! write (*,*) 'OUT OF TABLE: ',trim(tab%tableName),', i>imax-1, i / imax :', iindex,'/',tab%imaxTab, ii,jj,kk
      error = 1
    endif

    if (iindex < 2) then
      ! write (*,*) 'OUT OF TABLE: ',trim(tab%tableName),', i<2, i / imax :', iindex,'/',tab%imaxTab, ii,jj,kk
      error = 2
    endif


    jMin = 2
    jMax = tab%jmaxTab - 1

    if (tab%diagonalTable == 1)   then 
      jMin = 4
      jMax = tab%jmaxTabDiagonal - 3
    endif

    if (jindex > jMax) then
      ! write (*,*) 'OUT OF TABLE: ',trim(tab%tableName),', j>jmax, j / jmax :', jindex,'/',jMax, ii,jj,kk
      error = 3
    endif
    if (jindex < jMin) then
      ! write (*,*) 'OUT OF TABLE: ',trim(tab%tableName),', j<jmin, j / jmin :', jindex,'/',jMin, ii,jj,kk
      error = 4
    endif

    if (error > 0) return 

!   Check if inside the dome (tweo-phase)
    if (tab%diagonalTable == 1) then
      if (tab%val3(iindex+2,jindex-3) == 0.0_mytype) then 
        ! write(*,*) 'ERROR: ', tab%tableName, ' inside two phase'
        error = 5
      endif 
    endif

  end subroutine



  subroutine tableLookUp(error,ii,jj,kk, tab, dim1,dim2, & 
            val1,val2,val3,val4,val5,val6,val7,val8,val9,val10)

    implicit none
    real(mytype) :: dim1,dim2
    real(mytype), optional :: val1,val2,val3,val4,val5,val6,val7,val8,val9,val10
    integer :: i,j, jTable2d, factDiagTable
    integer :: ii,jj,kk
    real(mytype) :: error
    real(mytype) :: x0,x1,x2,x3
    real(mytype), dimension(4) ::fac_dim1, fac_dim2
    type (TABLE2D) :: tab

    i = int((dim1 - tab%dim1(1))/(tab%dim1(2)-tab%dim1(1))) + 1
    j = int((dim2 - tab%dim2(1))/(tab%dim2(2)-tab%dim2(1))) + 1

!     Check if the table is diagonal. If it is, then change the j index 
    jTable2d = j 
    if (tab%diagonalTable == 1)     jTable2d = j-i + 1

    call checkOutsideTable(i, jTable2d, tab, dim1, dim2, error, ii,jj,kk)

    if (error > 0) return 

    x0  = tab%dim1(i-1)
    x1  = tab%dim1(i  )
    x2  = tab%dim1(i+1)
    x3  = tab%dim1(i+2)

    fac_dim1(1) = (dim1-x1)*(dim1-x2)*(dim1-x3)/((x0-x1)*(x0-x2)*(x0-x3))
    fac_dim1(2) = (dim1-x0)*(dim1-x2)*(dim1-x3)/((x1-x0)*(x1-x2)*(x1-x3))
    fac_dim1(3) = (dim1-x0)*(dim1-x1)*(dim1-x3)/((x2-x0)*(x2-x1)*(x2-x3))
    fac_dim1(4) = (dim1-x0)*(dim1-x1)*(dim1-x2)/((x3-x0)*(x3-x1)*(x3-x2))

    x0 = tab%dim2(j-1)
    x1 = tab%dim2(j  )
    x2 = tab%dim2(j+1)
    x3 = tab%dim2(j+2)

    fac_dim2(1) = (dim2-x1)*(dim2-x2)*(dim2-x3)/((x0-x1)*(x0-x2)*(x0-x3))
    fac_dim2(2) = (dim2-x0)*(dim2-x2)*(dim2-x3)/((x1-x0)*(x1-x2)*(x1-x3))
    fac_dim2(3) = (dim2-x0)*(dim2-x1)*(dim2-x3)/((x2-x0)*(x2-x1)*(x2-x3))
    fac_dim2(4) = (dim2-x0)*(dim2-x1)*(dim2-x2)/((x3-x0)*(x3-x1)*(x3-x2))

    if(present(val1))  call lagrangeInterpol(val1,  tab%val1,  fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)
    if(present(val2))  call lagrangeInterpol(val2,  tab%val2,  fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)
    if(present(val3))  call lagrangeInterpol(val3,  tab%val3,  fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)
    if(present(val4))  call lagrangeInterpol(val4,  tab%val4,  fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)
    if(present(val5))  call lagrangeInterpol(val5,  tab%val5,  fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)
    if(present(val6))  call lagrangeInterpol(val6,  tab%val6,  fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)
    if(present(val7))  call lagrangeInterpol(val7,  tab%val7,  fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)
    if(present(val8))  call lagrangeInterpol(val8,  tab%val8,  fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)
    if(present(val9))  call lagrangeInterpol(val9,  tab%val9,  fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)
    if(present(val10)) call lagrangeInterpol(val10, tab%val10, fac_dim1, fac_dim2, i, jTable2d, tab%diagonalTable)

  end subroutine




  subroutine tableLookUp_Inverse(error, tab, dim1, tabvar, var, dim2)

    implicit none
    real(mytype) :: dim1, var, dim2
    integer :: i,j, jmax
    integer :: l, k
    real(mytype) :: error
    real(mytype) :: x0,x1,x2,x3
    real(mytype), dimension(4) ::fact, interpol
    type (TABLE2D) :: tab
    real(mytypeTable), dimension(:,:) :: tabvar

    i = int((dim1 - tab%dim1(1))/(tab%dim1(2)-tab%dim1(1))) + 1

    x0  = tab%dim1(i-1)
    x1  = tab%dim1(i  )
    x2  = tab%dim1(i+1)
    x3  = tab%dim1(i+2)

    fact(1) = (dim1-x1)*(dim1-x2)*(dim1-x3)/((x0-x1)*(x0-x2)*(x0-x3))
    fact(2) = (dim1-x0)*(dim1-x2)*(dim1-x3)/((x1-x0)*(x1-x2)*(x1-x3))
    fact(3) = (dim1-x0)*(dim1-x1)*(dim1-x3)/((x2-x0)*(x2-x1)*(x2-x3))
    fact(4) = (dim1-x0)*(dim1-x1)*(dim1-x2)/((x3-x0)*(x3-x1)*(x3-x2))

    jmax = tab%jmaxTabDiagonal

    j = 4
    do while (j<jmax-3)
      interpol = 0.0_mytype
      do k=-1,2
        do l=-1,2
          ! write(*,*) "index", i,j, k, l, i+l, j + k - tab%diagonalTable*l, tabvar(i+l, j + k - tab%diagonalTable*l)
          interpol(k+2) = interpol(k+2) + tabvar(i+l, j + k - tab%diagonalTable*l)*fact(l+2)
        enddo
      enddo
      ! write(*,*) "index", tab%imaxTab, jmax, i,j, var, interpol(2), interpol(3)
      if ((var-interpol(2))*(var-interpol(3)) .le. 0) exit
      j = j + 1
    enddo

    ! write(*,*) "index", tab%imaxTab, jmax, i,j, var, interpol(2), interpol(3)

    x0  = interpol(1)
    x1  = interpol(2)
    x2  = interpol(3)
    x3  = interpol(4)

    fact(1) = (var-x1)*(var-x2)*(var-x3)/((x0-x1)*(x0-x2)*(x0-x3))
    fact(2) = (var-x0)*(var-x2)*(var-x3)/((x1-x0)*(x1-x2)*(x1-x3))
    fact(3) = (var-x0)*(var-x1)*(var-x3)/((x2-x0)*(x2-x1)*(x2-x3))
    fact(4) = (var-x0)*(var-x1)*(var-x2)/((x3-x0)*(x3-x1)*(x3-x2))


    if (tab%diagonalTable == 1)  j = j+i-1 

    dim2 = 0.0_mytype
    do k=1,4
      dim2 = dim2 + tab%dim2(k-2+j)*fact(k)
    enddo

  end subroutine





  subroutine lagrangeInterpol(result, tabvar, fac_idir, fac_jdir, i, j, factDiagIndex)
    implicit none 
    integer :: i,j,k,l,factDiagIndex
    real(mytype), dimension(:) :: fac_idir, fac_jdir
    real(mytypeTable), dimension(:,:) :: tabvar
    real(mytype), dimension(4) :: interpol
    real(mytype) :: result

    interpol = 0.0_mytype

    do k=-1,2
      do l=-1,2
        ! write(*,*) "index", i,j, k, l, i+l, j + k - factDiagIndex*l, tabvar(i+l, j + k - factDiagIndex*l)
        interpol(k+2) = interpol(k+2) + tabvar(i+l, j + k - factDiagIndex*l)*fac_idir(l+2)
      enddo
    enddo

    result = 0.0_mytype
    do k=1,4
      result = result + interpol(k)*fac_jdir(k)
    enddo

  end subroutine



end module



