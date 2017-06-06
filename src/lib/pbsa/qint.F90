#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine qint here]
subroutine qint(l,m,n,h,hx,hy,hz,xs,ys,zs,maxirr,nq,index,index2,cirreg, &
              wcoe,wycoe,wzcoe,q0p)

   use iim_util
   implicit none

   ! Passed variables

   integer l,m,n,maxirr,nq
   _REAL_  h,hx,hy,hz,xs,ys,zs
   integer index(l,m,n), index2(l,m,n)
   _REAL_  cirreg(15, maxirr)
   _REAL_  wcoe(maxirr,nq),wycoe(maxirr, nq),wzcoe(maxirr, nq)
   _REAL_  q0p(maxirr)

   !Local variables

   integer i,j,k,ir
   _REAL_  q,qy,qz

   do i=1, l
   do j=1, m
   do k=1, n
      if (index(i,j,k) > 1 .and. index(i,j,k) < 5) then
         ir = index2(i, j, k)
         call qint2(l,m,n,h,hx,hy,hz,xs,ys,zs,ir,maxirr,nq,index2,cirreg,qp,&
            wcoe,wycoe,wzcoe,q,qy,qz)
         q0p(ir) = q
         qyp(ir) = qy
         qzp(ir) = qz
      end if
   end do
   end do
   end do

end subroutine qint
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine qint2 here]
subroutine qint2(l,m,n,h,hx,hy,hz,xs,ys,zs,ir, maxirr,nq,index2,cirreg,qp,&
              wcoe,wycoe,wzcoe,q,qy,qz)

   implicit none
   integer,parameter :: nqq = 27

   ! Passed variables

   integer l,m,n,ir,maxirr,nq
   _REAL_  xs,ys,zs,h,hx,hy,hz
   integer index2(l,m,n)
   _REAL_  cirreg(15, maxirr)
   _REAL_  wcoe(maxirr,nq),wycoe(maxirr, nq),wzcoe(maxirr, nq)
   _REAL_  qp(maxirr),t(3,3)
   _REAL_  q,qy,qz
   _REAL_  itemp(nqq)

   !Local variables

   integer if0,i0,j0,k0,nsub,i1,j1,k1,i2,j2,k2,idis,ii,ip,nn2
   _REAL_  hmax,alf,x1,y1,z1,xyy,xzz,xyz,x2,y2,z2,dis,x3,y3,z3

   hmax=h
   if (nq /= nqq) then
      write(*,*) "   Please change nqq in qint2() s.t. nqq=nq in main()!"
      stop
   end if

   alf = 10.1*hmax
   if0 = int(alf/hmax) + 1

   q  = 0.0
   qy = 0.0
   qz = 0.0

   x1     = cirreg( 1, ir)
   y1     = cirreg( 2, ir)
   z1     = cirreg( 3, ir)

   xyy    = cirreg( 4, ir)
   xzz    = cirreg( 5, ir)
   xyz    = cirreg( 6, ir)

   t(1,1) = cirreg( 7, ir)
   t(1,2) = cirreg( 8, ir)
   t(1,3) = cirreg( 9, ir)
   t(2,1) = cirreg(10, ir)
   t(2,2) = cirreg(11, ir)
   t(2,3) = cirreg(12, ir)
   t(3,1) = cirreg(13, ir)
   t(3,2) = cirreg(14, ir)
   t(3,3) = cirreg(15, ir)

   i0 = nint((x1-xs)/hx)
   j0 = nint((y1-ys)/hy)
   k0 = nint((z1-zs)/hz)

   nsub = 0
   do i2 = 0, if0          !# Starting from the center
      do i1 = i0-i2, i0+i2
      do j1 = j0-i2, j0+i2
      do k1 = k0-i2, k0+i2
         if (i1 < 1 .or. j1 < 1 .or. k1 < 1 .or. i1 > l .or. j1 > m .or. k1 > n) goto 555

         idis = abs(i1-i0) + abs(j1-j0) + abs(k1-k0)
         nn2 = index2(i1,j1,k1)

         if (idis == i2 .and. nn2 > 0) then
            x2 = cirreg(1, nn2)
            y2 = cirreg(2, nn2)
            z2 = cirreg(3, nn2)

            do ii = 1, nsub
               ip = itemp(ii)
               x3 = cirreg(1, ip)
               y3 = cirreg(2, ip)
               z3 = cirreg(3, ip)
               call distan(x2,y2,z2,x3,y3,z3, dis)
               if (dis < 0.1*hmax) then
                  goto 555
               end if
            end do

            call distan(x1,y1,z1,x2,y2,z2, dis)
            if ( dis < alf ) then
               nsub = nsub + 1
               itemp(nsub) = nn2
               ip = index2(i1,j1,k1)
               q  = q + wcoe(ir, nsub)*qp(ip)
               qy = qy + wycoe(ir, nsub)*qp(ip)
               qz = qz + wzcoe(ir, nsub)*qp(ip)
               if ( nsub .eq. nq) goto 577
            end if
         end if

         555 continue
      end do  ! k1 = k0-i2, k0+i2
      end do  ! j1 = j0-i2, j0+i2
      end do  ! i1 = i0-i2, i0+i2
   end do  ! i2=0, if0
   577 continue

end subroutine qint2
