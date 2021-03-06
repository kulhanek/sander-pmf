#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine coed20 here]
subroutine coed20(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,maxirr,nq,index2,cirreg, &
              wcoe, wxcoe, wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)

   use iim_util
   implicit none

   ! Passed variables

   integer l,m,n,nirreg,maxirr,nq
   integer index2(l,m,n)
   _REAL_  xs,ys,zs,h,hx,hy,hz
   _REAL_  cirreg(15, maxirr)
   _REAL_  wcoe(maxirr,nq),wxcoe(maxirr, nq),wycoe(maxirr, nq)
   _REAL_  wzcoe(maxirr,nq),wxxcoe(maxirr, nq),wyycoe(maxirr, nq)
   _REAL_  wzzcoe(maxirr, nq),wxycoe(maxirr, nq)
   _REAL_  wxzcoe(maxirr, nq),wyzcoe(maxirr, nq)

   ! Local variables

   integer,parameter :: nqq = 27
   integer i0,j0,k0,ir,isvd,i,j,k,i1,j1,k1,ip,job,ii,&
           ij,jj,nsub,nn2,i2,if0,idis,inf
   _REAL_  alf,x1,y1,z1,x2,y2,z2,x3,y3,z3,xyy,xzz,xyz,dis,hmax
   _REAL_  w1(10, nqq), w4(10)
   _REAL_  sd(11), uw(10,10), v(nqq, nqq)
   _REAL_  ew(nqq)
   _REAL_  ew1(nqq), ew2(nqq), ew3(nqq)
   _REAL_  ew4(nqq), ew5(nqq), ew6(nqq)
   _REAL_  ew7(nqq), ew8(nqq), ew9(nqq), ew10(nqq)
   _REAL_  w31(nqq), w32(nqq), w33(nqq)
   _REAL_  w34(nqq), w35(nqq), w36(nqq)
   _REAL_  w37(nqq), w38(nqq), w39(nqq), w310(nqq)
   _REAL_  sd2(10), work(1000)
   _REAL_  t(3,3)
   _REAL_  tempx(3),tempy(3)
   integer itemp(nqq)

   hmax=h
   isvd = 1

   if (nq /= nqq) then
      write(*,*) "   Please change nqq in coedef() s.t. nqq=nq in main()!"
      stop
   end if

   alf = 18.1*hmax
   if0 = int(alf/hmax) + 1

   do i=1, nirreg
      do j=1, nq
         wcoe(i,j)   = 0.0
         wxcoe(i,j)  = 0.0
         wycoe(i,j)  = 0.0
         wzcoe(i,j)  = 0.0
         wxxcoe(i,j) = 0.0
         wyycoe(i,j) = 0.0
         wzzcoe(i,j) = 0.0
         wxycoe(i,j) = 0.0
         wxzcoe(i,j) = 0.0
         wyzcoe(i,j) = 0.0
      end do
   end do

   ! ----- for each projection point

   do ir = 1, nirreg

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

      ! -------- initialize

      do jj=1, nq
         do ii=1, nq
            v(ii, jj) = 0.0      !# the V matrix in SVD: U*S*V
         end do
         do ii=1, 20
            w1(ii, jj) = 0.0     !# the coefficient matrix to be SVDed
         end do
      end do

      ! -------- Generate the matrix

      nsub = 0
      do i2=0, if0          !# Starting from the center
         do i1 = i0-i2, i0+i2
         do j1 = j0-i2, j0+i2
         do k1 = k0-i2, k0+i2
            if (i1 < 1 .or. j1 < 1 .or. k1 < 1 .or. i1 > l .or. j1 > m .or. k1 > n) goto 55

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
                     goto 55
                  end if
               end do

               call distan(x1,y1,z1,x2,y2,z2, dis)
               if ( dis < alf ) then
                  nsub = nsub + 1
                  itemp(nsub) = nn2
                  tempx(1) = x2 - x1
                  tempx(2) = y2 - y1
                  tempx(3) = z2 - z1
                  call matvec(3,3,t, tempx, tempy)
                  w1(1, nsub)  = 1.0
                  w1(2, nsub)  = tempy(1)
                  w1(3, nsub)  = tempy(2)
                  w1(4, nsub)  = tempy(3)
                  w1(5, nsub)  = 0.5*tempy(1) * tempy(1)
                  w1(6, nsub)  = 0.5*tempy(2) * tempy(2)
                  w1(7, nsub)  = 0.5*tempy(3) * tempy(3)
                  w1(8, nsub)  = tempy(1) * tempy(2)
                  w1(9, nsub)  = tempy(1) * tempy(3)
                  w1(10, nsub) = tempy(2) * tempy(3)
                  if ( nsub .eq. nq) goto 77
               end if
            end if  ! (idis == i2 .and. nn2 > 0)
            55 continue
         end do  ! k1 = k0-i2, k0+i2
         end do  ! j1 = j0-i2, j0+i2
         end do ! i1 = i0-i2, i0+i2
      end do ! i2=0, if0
      77 continue

      if (nsub < 12) then
         write(*,*) "   nsub = ",nsub,"alf is too small in coede20f()!"
         stop
      end if

      ! -------- Call least square routine

      if (isvd == 1) then
         job = 11
         call dsvdc(w1,10,10,nsub,sd,ew,uw,10,v,nq,w4,job,inf)
      else
         call dgesvd('A','A',10,nsub,w1,10,sd2,uw,10,v,nq, &
               work, 1000, inf)

         do ij = 1, 10
            sd(ij) = sd2(ij)
         end do
      end if

      if (inf /= 0) then
         write(*,*) inf, " - ssvdc() or sgesvd() failed in coedef!"
         stop
      end if

      do i1=1, 10
         if (abs(sd(i1)) > 1.0e-14 ) then
            ew1(i1)=uw(1, i1)/sd(i1)
            ew2(i1)=uw(2, i1)/sd(i1)
            ew3(i1)=uw(3, i1)/sd(i1)
            ew4(i1)=uw(4, i1)/sd(i1)
            ew5(i1)=uw(5, i1)/sd(i1)
            ew6(i1)=uw(6, i1)/sd(i1)
            ew7(i1)=uw(7, i1)/sd(i1)
            ew8(i1)=uw(8, i1)/sd(i1)
            ew9(i1)=uw(9, i1)/sd(i1)
            ew10(i1)=uw(10, i1)/sd(i1)
         else
            ew1(i1)= 0.0
            ew2(i1)= 0.0
            ew3(i1)= 0.0
            ew4(i1)= 0.0
            ew5(i1)= 0.0
            ew6(i1)= 0.0
            ew7(i1)= 0.0
            ew8(i1)= 0.0
            ew9(i1)= 0.0
            ew10(i1)= 0.0
         end if
      end do

      ! -------- w31(i),w32(i),......,w36(i) are the solutions

      do i1=1, nsub
         w31(i1) = 0.0
         w32(i1) = 0.0
         w33(i1) = 0.0
         w34(i1) = 0.0
         w35(i1) = 0.0
         w36(i1) = 0.0
         w37(i1) = 0.0
         w38(i1) = 0.0
         w39(i1) = 0.0
         w310(i1) = 0.0

         if (isvd == 1) then
            do j1=1,10
               w31(i1) = w31(i1) + v(i1,j1)*ew1(j1)
               w32(i1) = w32(i1) + v(i1,j1)*ew2(j1)
               w33(i1) = w33(i1) + v(i1,j1)*ew3(j1)
               w34(i1) = w34(i1) + v(i1,j1)*ew4(j1)
               w35(i1) = w35(i1) + v(i1,j1)*ew5(j1)
               w36(i1) = w36(i1) + v(i1,j1)*ew6(j1)
               w37(i1) = w37(i1) + v(i1,j1)*ew7(j1)
               w38(i1) = w38(i1) + v(i1,j1)*ew8(j1)
               w39(i1) = w39(i1) + v(i1,j1)*ew9(j1)
               w310(i1) = w310(i1) + v(i1,j1)*ew10(j1)
            end do
         else
            do j1=1,10
               w31(i1) = w31(i1) + v(j1,i1)*ew1(j1)
               w32(i1) = w32(i1) + v(j1,i1)*ew2(j1)
               w33(i1) = w33(i1) + v(j1,i1)*ew3(j1)
               w34(i1) = w34(i1) + v(j1,i1)*ew4(j1)
               w35(i1) = w35(i1) + v(j1,i1)*ew5(j1)
               w36(i1) = w36(i1) + v(j1,i1)*ew6(j1)
               w37(i1) = w37(i1) + v(j1,i1)*ew7(j1)
               w38(i1) = w38(i1) + v(j1,i1)*ew8(j1)
               w39(i1) = w39(i1) + v(j1,i1)*ew9(j1)
               w310(i1) = w310(i1) + v(j1,i1)*ew10(j1)
            end do
         end if

         wcoe(ir,i1)   = w31(i1)
         wxcoe(ir,i1)  = w32(i1)
         wycoe(ir,i1)  = w33(i1)
         wzcoe(ir,i1)  = w34(i1)
         wxxcoe(ir,i1) = w35(i1)
         wyycoe(ir,i1) = w36(i1)
         wzzcoe(ir,i1) = w37(i1)
         wxycoe(ir,i1) = w38(i1)
         wxzcoe(ir,i1) = w39(i1)
         wyzcoe(ir,i1) = w310(i1)
      end do  !  i1=1, nsub

   end do  ! ir = 1, nirreg

end subroutine coed20
