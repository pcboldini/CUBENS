! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! standard units module

module io_std_units
   use decomp_2d
   implicit none
   ! specify the output mode
#ifdef OUTPUT_UNIT
   ! standard output stream
   integer(mytype), parameter :: stdout = 6                                     
#else
   ! no screen output but in a separate file
   integer(mytype), parameter :: stdout = 7                                     
#endif   
   
end module io_std_units
