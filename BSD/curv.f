#define _REAL_ double precision

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine curv here]
subroutine curv(x1,y1,z1,i0,j0,k0,x,y,z,phi,t, xyy,xzz,xyz)
   implicit _REAL_ (a-h,o-z)

   !       ----------------------------------------------
   !       Compute curvature at the projection (x1,y1,z1)
   !       ----------------------------------------------

   common /lmn/l, m, n, nirreg
   common /hxyz/hx,hy,hz,hmax

   dimension x(0:l+1), y(0:m+1), z(0:n+1)
   dimension phi(0:l+1,0:m+1,0:n+1)

   dimension t(3,3)

   call grtopr(x,y,z,x1,y1,z1,i0,j0,k0,phi,1,1,1, &
         ph0,phx,phy,phz,phxx,phyy,phzz,phxy,phxz,phyz)

   phieps = t(1,1)*phx+t(1,2)*phy+t(1,3)*phz

   if ( abs(phieps) < 1.0e-20) then
      write(*,*) "   Error: phieps = 0.0 in curv()!"
      stop
   end if

   tmp1 = t(2,1)*phxx+t(2,2)*phxy+t(2,3)*phxz
   tmp2 = t(2,1)*phxy+t(2,2)*phyy+t(2,3)*phyz
   tmp3 = t(2,1)*phxz+t(2,2)*phyz+t(2,3)*phzz

   tmp = t(2,1)*tmp1+t(2,2)*tmp2+t(2,3)*tmp3
   xyy = -tmp/phieps

   tmp1 = t(3,1)*phxx+t(3,2)*phxy+t(3,3)*phxz
   tmp2 = t(3,1)*phxy+t(3,2)*phyy+t(3,3)*phyz
   tmp3 = t(3,1)*phxz+t(3,2)*phyz+t(3,3)*phzz

   tmp = t(3,1)*tmp1+t(3,2)*tmp2+t(3,3)*tmp3
   xzz = -tmp/phieps

   tmp = t(2,1)*tmp1+t(2,2)*tmp2+t(2,3)*tmp3
   xyz = -tmp/phieps

   return
end subroutine curv 

