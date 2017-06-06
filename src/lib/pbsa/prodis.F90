#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine prodis here]
subroutine prodis(l,m,n,nirreg,bcopt,bi,bo,atmfirst,atmlast,maxirr,cirreg)

   use iim_util
   implicit none

   ! passed variables

   integer l, m, n, nirreg
   integer bcopt, atmfirst, atmlast, maxirr
   _REAL_ cirreg(15, maxirr)
   _REAL_ t(3,3)
   _REAL_ bi,bo

   ! local variables

   integer ir
   _REAL_ xx, yy, zz

   ! external functions

   _REAL_ fw, fq

   do ir = 1, nirreg
      xx     = cirreg( 1, ir)
      yy     = cirreg( 2, ir)
      zz     = cirreg( 3, ir)
      t(1,1) = cirreg( 7, ir)
      t(1,2) = cirreg( 8, ir)
      t(1,3) = cirreg( 9, ir)
      t(2,1) = cirreg(10, ir)
      t(2,2) = cirreg(11, ir)
      t(2,3) = cirreg(12, ir)
      t(3,1) = cirreg(13, ir)
      t(3,2) = cirreg(14, ir)
      t(3,3) = cirreg(15, ir)
      wp(ir) = fw(bcopt,bi,bo,atmfirst,atmlast,xx,yy,zz)
      qp(ir) = fq(bcopt,bi,bo,atmfirst,atmlast,xx,yy,zz,t)
   end do

end subroutine prodis
