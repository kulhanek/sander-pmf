#include "copyright.h"
#include "../include/dprec.fh"

module iim_util

   _REAL_, allocatable :: wp(:),qp(:)
   _REAL_, allocatable :: qyp(:),qzp(:)
   _REAL_, allocatable :: wyp(:),wzp(:)
   _REAL_, allocatable :: wyyp(:),wzzp(:),wyzp(:)

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ local copy of matvec
subroutine matvec(m,n,a,x,y)

   implicit none
   integer m,n,i,j
   _REAL_  a(m,n), x(n),y(m)

   do i=1,m
      y(i) = 0.0
      do j=1,n
         y(i) = y(i) + a(i,j)*x(j)
      end do
   end do

end subroutine matvec

end module iim_util
