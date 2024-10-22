! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
! math routines module

module mod_math
   use decomp_2d
   implicit none
   ! definition of double precision PI
   real(mytype), parameter :: pi_const    = 4.0_mytype*atan(1.0_mytype)


contains
   ! calculation of cubic roots for EoS
   subroutine cubic_root(A,B,C,D,x)
      !$acc routine seq
      use decomp_2d, only: mytype
      implicit none 
      real(mytype), intent(in)  :: A,B,C,D
      real(mytype), intent(OUT) :: x
      real(mytype) :: thrd,secd,prec
      real(mytype) :: delta_0, delta_1, disc
      real(mytype) :: x1_atan, x2_atan, x3_atan
      complex*16 :: comp_fact, c_fact, x1, x2, x3
      thrd = 1.0_mytype/3.0_mytype
      secd = 1.0_mytype/2.0_mytype
      prec = 1.0e-4_mytype
      comp_fact=secd*cmplx(-1,3.0_mytype**secd)
      delta_0=B**2-3*A*C
      delta_1=2*B**3-9*A*B*C+27*A**2*D
      disc=( delta_1**2 -4*delta_0**3 )**secd
      c_fact=cmplx((disc+delta_1)/2,0.0_mytype)**thrd
      x1=-thrd/A*(B+c_fact+delta_0/c_fact)
      x2=-thrd/A*(B+c_fact*comp_fact+delta_0/c_fact/comp_fact)
      x3=-thrd/A*(B+c_fact*comp_fact**2+delta_0/c_fact/comp_fact**2)
      x1_atan=atan2(aimag(x1),real(x1))
      x2_atan=atan2(aimag(x2),real(x2))
      x3_atan=atan2(aimag(x3),real(x3))
      if  ( (abs(x1_atan) .LE. prec) .AND. (abs(aimag(x1)) .LE. prec) ) then
         x=real(x1)
      elseif  ( (abs(x2_atan) .LE. prec) .AND. (abs(aimag(x2)) .LE. prec) ) then
         x=real(x2)
      elseif  ( (abs(x3_atan) .LE. prec) .AND. (abs(aimag(x3)) .LE. prec) ) then
         x=real(x3)
      endif
   end subroutine


! spline function for mesh cubic interpolation: calculation of the second derivative
   subroutine spline(x, y, n, y2)
      !$acc routine seq
      use decomp_2d
      implicit none
      integer   i, k, n, nmax
      real(mytype)    yp1, ypn, x(n), y(n), y2(n), p(n), p2(n), qn, sig(n), un, u(n)
      y2(1) = 0.0_mytype
      u(1)  = 0.0_mytype
      p(1)  = 1.0_mytype 
      sig(1) = 0.0_mytype
      do i=2, n-1
         sig(i)=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p(i)=sig(i)*y2(i-1)+2
         y2(i)=(1.0_mytype*sig(i)-1.0_mytype)/p(i)
         u(i)=(6.0_mytype*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) / & 
              (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig(i)*u(i-1))
      enddo
      u(n) = u(n-1)
      p(n) = p(n-1)
      u=u/p
      qn=0.0_mytype
      un=0.0_mytype
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1)
      do k=n-1, 1, -1
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
   end subroutine


! spline function for mesh cubic interpolation: actual interpolation with given mesh second derivatives
subroutine splint(xa,ya,y2a,n,x,y)
   !$acc routine seq
   use decomp_2d
   implicit none
   integer n,k,khi,klo
   real(mytype) x,y,xa(n),y2a(n),ya(n), a,b,h
   klo=1
   khi=n
1  if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
         khi=k
      else
         klo=k
      endif
      goto 1
   endif
   h = xa(khi)-xa(klo)
   if (h.eq.0.) stop 'bad xa input in splint'
   a = (xa(khi)-x)/h
   b = (x-xa(klo))/h
   y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6
   return
end subroutine
! numerical integration: trapezoidal method
subroutine trapzint(n,xa,ya,y)
   use decomp_2d
   implicit none
   integer n
   real(mytype) x,y,xa(n),ya(n)
   y=sum((ya(1+1:n-0) + ya(1+0:n-1))*(xa(1+1:n-0) - xa(1+0:n-1)))/2
   return
end subroutine


end module
