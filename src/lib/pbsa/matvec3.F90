#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Caculate A*X with the Un inputted,which is [beta*un] + bf  ]
subroutine matvec3(xm,ym,zm,h,nbnd,nind3,x,y,z,xs,ys,zs,phi,u,f,&
             cirreg,&
             unj,unjf,fvec,bf,epsin,epsout,index,index2,&
             q0p,w0p,&
             wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe,&
             wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2,&
             accept,bv,icall,solvopt,norm,inorm )

   use iim_util
   use aug_solver   ! change interface names, 'personalized' solver used.

   ! input unj, output fvec
   ! nind3 is the  reduced number of  all the irregular points
   ! note unjf and  bf are passed in u&u0 are not
   ! note this is two sides aug approaches.

   implicit none
   integer, parameter :: nq=27
   logical alive
   integer :: status = 0

   ! passed variables

   ! nbnd,nind3              the initial dimension; reduced dimension of aug
   ! x(),y(),z()             coordinates of grid points
   ! xs,ys,zs                staring crds of the box
   ! phi,u                   lvlset function; potential of the whole grids
   ! cirreg                  geometrical parameters at irregular points
   ! index()                 flag of all grid points with 1-5
   ! index2()                number index for irregular grid points
   ! wp,qp                   [u]and[betaUn]for the interface
   ! unj(),unjf()            [Un] for irregular points, unjf is the input
   ! bf(), fvec()            rhs of Ag=b; output of matrix-vector multiply
   ! w0p,wyp,wzp,w0p,        tangential derivatives of jump conditions
   ! wyp,wzp
   ! wcoe,wxcoe,wycoe,wzcoe  coefficients obtained by SVD
   ! wxxcoe,wyycoe wzzcoe
   ! wxycoe,wxzcoe,wyzcoe
   ! bv                      initial rhs of the linear eqns of the whole grids
   ! solvopt,icall           solvopt for linear solvers above, debuging variable
   ! c,c2                    linear coefficients
   ! accept                  convergence criteria of the linear solvers
   ! sss1,sss2               dummy arrays

   ! local variables

   ! xp,yp,zp                projection crds on interface of the irregular pnts
   ! eps[x-z]                the coefs of grid pnts, which is 1 for AUG method
   ! f(i,j,k)                rhs of the eqns for the whole grid for AUG method
   ! xyy,xzz,xyz             derivatives of the geometry on projection pnts
   ! t                       trans matrix between local crds and grid crds
   ! unin                    derivative of u from inner side on proj pnts

   integer xm,ym,zm,nbnd,nind3,n3,n1,icall
   integer index(1:xm,1:ym,1:zm),index2(1:xm,1:ym,1:zm)
   _REAL_ xyz,xyy,xzz,gox,goy,goz,h,accept,epsout,epsin
   _REAL_  q0p(nbnd),w0p(nbnd)

   _REAL_  wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd,nq),wzcoe(nbnd,nq),&
           wxxcoe(nbnd, nq),wyycoe(nbnd, nq),wzzcoe(nbnd, nq),&
           wxycoe(nbnd, nq),wxzcoe(nbnd, nq),wyzcoe(nbnd,nq)
   _REAL_  sss1(nbnd),sss2(nbnd)
   _REAL_  c(xm,ym,zm,7),c2(nbnd, 27)
   _REAL_  t(3,3)
   _REAL_  bv(xm,ym,zm)
   _REAL_  x(0:xm+1),y(0:ym+1),z(0:zm+1)
   _REAL_  phi(0:xm+1,0:ym+1,0:zm+1)
   _REAL_  u(1:xm,1:ym,1:zm),f(xm,ym,zm)
   _REAL_  unj(1:nbnd),cirreg(1:15,1:nbnd),&
           unjf(1:nbnd),fvec(1:nbnd),bf(1:nbnd)

   ! variables to use mg

   _REAL_  inorm,norm,dummy
   integer itn , solvopt
   _REAL_  epsx,epsy,epsz
   _REAL_  iv(1:xm*ym*zm)
   _REAL_  xso(xm*ym*zm+2*xm*ym)

   _REAL_ xs,ys,zs,xf,yf,zf,hx,hy,hz,hmax,bi,bo
   _REAL_ uu,dudx,dudy,dudz,unin,unout
   integer l, m, n, nirreg

   !sgni,j,k is not used for now... XL

   integer :: sgni=0
   integer :: sgnj=0
   integer :: sgnk=0
   integer i0,j0,k0,nn1

   !WMBS - for use in neutralizing charge on f for periodic solvers

   integer cntirreg !count of the number of irregular grid nodes
   _REAL_ fsum !holds the average charge per irregular grid node to
               !be neutralized

   integer i,j,k,IFAIL,IFAIL2,ii,ir,nc,nz_num
   _REAL_ beta_max,rhs
   _REAL_ coe2(27)
   _REAL_ fw

   _REAL_ xp,yp,zp

   _REAL_ qpsav(nbnd)

   l = xm; m = ym; n = zm; nirreg = nind3
   hx=h;hy=h;hz=h

   epsx=1.0d0;epsy=1.0d0;epsz=1.0d0
   bi=epsin;bo=epsout

   cntirreg = 0

   do i = 1, nind3
      unj(i) = unjf(i)
      qpsav(i)=qp(i)
      qp(i)=unjf(i)
   end do

   ! calculate the first and second derivatives of the jump conditions in the
   ! surface tangential directions in the local coordinate system
   ! step 1: this is the first derivatives of w g for  irregular points only

   call coed20(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
           wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)

   call qint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg, &
           wcoe,wycoe,wzcoe,q0p)
!  call qint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,unj, &
!            wcoe,wycoe,wzcoe,q0p,qyp,qzp)

   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
           wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,wyp,wzp,wyyp,wzzp,wyzp)


   ! step 2: this is the second derivatives ( when secodary order, should be used )

   call coed6(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
           wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)

   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
           wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,sss1,sss2,wyyp,wzzp,wyzp)

   ! setting up linear system coefficient matrix
   ! the current version uses irregular points on both sides

   beta_max=1.0
   nz_num = 0

   do k = 1, n
   do j = 1, m
   do i = 1, l

      ! inside regular points

      if ( index(i,j,k) == 1 ) then

         ! set up the 7-band laplassian operator in uniform beta_max

         do ii = 2, 7
            c(i,j,k,ii) =1.0d0/h/h             !epsin/h/h
         end do
         if ( i == 1 ) c(i,j,k,2) = 0.d0
         if ( i == l ) c(i,j,k,3) = 0.d0
         if ( j == 1 ) c(i,j,k,4) = 0.d0
         if ( j == m ) c(i,j,k,5) = 0.d0
         if ( k == 1 ) c(i,j,k,6) = 0.d0
         if ( k == n ) c(i,j,k,7) = 0.d0
         c(i,j,k,1) = 6.0d0/h/h ! epsin/h/h  XP: this is still for general IIM, epsin /= epsout
         f(i,j,k) = bv(i,j,k)/h/h/h/epsin*FOURPI
         do ii = 1 , 7
            if ( abs(c(i,j,k,ii)) > 1.d-10 ) nz_num = nz_num + 1
         end do

      ! outside regular points

      else if (index(i,j,k) == 5 ) then

         do ii = 2, 7
            c(i,j,k,ii) = 1.0d0/h/h          !epsout/h/h
         end do
         if ( i == 1 ) c(i,j,k,2) = 0.d0
         if ( i == l ) c(i,j,k,3) = 0.d0
         if ( j == 1 ) c(i,j,k,4) = 0.d0
         if ( j == m ) c(i,j,k,5) = 0.d0
         if ( k == 1 ) c(i,j,k,6) = 0.d0
         if ( k == n ) c(i,j,k,7) = 0.d0
         c(i,j,k,1) = 6.0d0/h/h ! epsout/h/h ! XP: this is still for general IIM, epsin /= epsout
         f(i,j,k) = bv(i,j,k)/h/h/h/epsout*FOURPI
         do ii = 1, 7
            if ( c(i,j,k,ii ) > 1.d-10 ) nz_num=nz_num+1
         end do

      ! irregular points

      else
         cntirreg = cntirreg+1  !WS
         coe2=0.0d0

         ! XP: Note the interface has changed. b_in, b_out, and bi bo are all passed in.

         bo=bi
         call irre31(l,m,n,h,hx,hy,hz,IFAIL, &
                 i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
                 nq,nbnd,index2,cirreg,coe2,rhs)
!        call irre31(l,m,n,h,hx,hy,hz,IFAIL, &
!                i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
!                qyp,qzp,wyp,wzp,wyyp,wzzp,wyzp,  &
!                nq,nbnd,index2,cirreg,wp,unj,coe2,rhs)

         if (IFAIL.gt.10) call irre32(l,m,n,h,hx,hy,hz,IFAIL2, &
                 i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
                 q0p,w0p, &
                 nq,nbnd,index2,cirreg,coe2,rhs)
!        if (IFAIL.gt.10) call irre32(l,m,n,h,hx,hy,hz,IFAIL2, &
!                i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
!                q0p,qyp,qzp,w0p,wyp,wzp,wyyp,wzzp,wyzp, &
!                nq,nbnd,index2,cirreg,wp,unj,coe2,rhs)

         ir = index2(i,j,k)
         bo=epsout
         f(i,j,k) = -rhs            !*(-6.0d0)/coe2(14)

         ! Out of the 27 neighbors, 7 is nonzero and contribute

         do nc = 1, 27
            c2(ir,nc) = coe2(nc)
         end do
         do ii = 1, 27
            if ( abs(c2(ir,ii)) > 1.d-10 ) nz_num = nz_num + 1
         end do
      end if
   end do
   end do
   end do

! XP: Outputs check for the setup of the linear eqns entering the linear system solver
!if(icall == 2 ) then
!   write(7001,*) l,m,n,nbnd,nz_num
!   write(7002,"(7f20.6)") c
!   write(7007,"(7f20.6)") c2
!   write(7003,*) index,index2
!   write(7004,*) f               !error
!   write(7006,"(7f20.6)") u
!   write(7005,*) u0
!   write(7005,*) accept,icall
!   write(7008,*) "h=",h,"epsin=",epsin,"epsout=",epsout
!end if
!if(icall == 1 ) then
!   write(6009,*) bv
!   write(6001,*) l,m,n,nbnd,nz_num
!   write(6002,"(7f20.6)") c
!   write(6007,"(7f20.6)") c2
!   write(6003,*) index,index2
!   write(6004,*) f               !error
!   write(6006,"(7f20.6)") u
!   write(6005,*) u0
!   write(6005,*) accept
!   write(6008,*) "h=",h,"epsin=",epsin,"epsout=",epsout
!end if

   f=f*h*h
   dummy=0.0d0;iv=0.0d0; !use iccg
   call init_param(l,m,n,l*m,l*m*n,10000,dummy,accept, 0.0d0,1.0d0,h,1.9d0)
   call allocate_array(solvopt)
   xso(:) = 0.d0
   call init_array(solvopt,epsx,epsy,epsz,f,iv,xso)

   if (solvopt /= 7 ) then !.and. solvopt/=5) then
      call pb_iccg(u,xso)
   !else if ( solvopt==7) then
   !   !WMBS - Need to neutralize any excess charge before calling fft.
   !   !FFT will do this internally, but it won't show up in f since it
   !   !will switch to an internal charge grid.
   !   if (abs(sum(f)/(xm*ym*zm)) > 1.0d-8) then
   !      fsum = sum(f)/(xm*ym*zm)
   !      do i=1,l; do j=1,m; do k=1,n
   !          f(i,j,k) = f(i,j,k) - fsum
   !      end do; end do; end do;
   !      flush(6)
   !   end if
   !   call pb_fftsolv(f,2,xm,ym,zm,u,h,64)
   !   u = u*h
   end if

   itn = l_itn
   inorm = l_inorm
   norm = l_norm

   call deallocate_array(solvopt)

   n3=10! second order interplation
   do k=1,n
   do j=1,m
   do i=1,l
      nn1=index2(i,j,k)

      if (nn1 .gt. 0) then

         ! projection crds and geometrical info

         xp     = cirreg( 1, nn1)
         yp     = cirreg( 2, nn1)
         zp     = cirreg( 3, nn1)

         xyy    = cirreg( 4, nn1)
         xzz    = cirreg( 5, nn1)
         xyz    = cirreg( 6, nn1)

         t(1,1) = cirreg( 7, nn1)
         t(1,2) = cirreg( 8, nn1)
         t(1,3) = cirreg( 9, nn1)
         t(2,1) = cirreg(10, nn1)
         t(2,2) = cirreg(11, nn1)
         t(2,3) = cirreg(12, nn1)
         t(3,1) = cirreg(13, nn1)
         t(3,2) = cirreg(14, nn1)
         t(3,3) = cirreg(15, nn1)

         !the nearest grid point

         i0 = nint((xp - x(0))/h)
         j0 = nint((yp - y(0))/h)
         k0 = nint((zp - z(0))/h)

         bo=bi
         call twosided(l,m,n,bi,bo,n3,i0,j0,k0,&
                 nbnd,xp,yp,zp,uu,dudx,dudy,dudz,&
                 xyy,xzz,xyz,t,u,phi,x(0),y(0),z(0),h,&
                 nn1)
!        call twosided(l,m,n,bi,bo,n3,i0,j0,k0,&
!                nbnd,xp,yp,zp,uu,dudx,dudy,dudz,&
!                wp,wyp,wzp,wyyp,wzzp,wyzp,unj,qyp,qzp, &
!                xyy,xzz,xyz,t,u,phi,x(0),y(0),z(0),h,&
!                nn1)
         bo=epsout

         ! call interp output Un^-, dudx in the local crds

         unin = dudx
         unout = unin + unj(nn1)

         ! Preconditioning, which results in smaller residual
         ! unout = (qp(nn1)-unj(nn1))/79.0d0 , as we are using 1:80

         fvec(nn1)=epsout*unout-epsin*unin-qpsav(nn1)
         fvec(nn1)=fvec(nn1)+bf(nn1)
      end if
   end do
   end do
   end do

   do i = 1, nind3
      qp(i) = qpsav(i)
   end do

end subroutine matvec3
