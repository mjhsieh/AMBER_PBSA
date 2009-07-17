! <compile=optimized>
#include "copyright.h"
#include "is_copyright.h"
#include "dprec.h"
#include "is_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for molecular dielectric map assignment.
subroutine pb_exmol( pbverbose,ifcap )

   use poisson_boltzmann    
   use solvent_accessibility

   implicit none

   ! Passed variables
 
   logical pbverbose
   integer ifcap
 
   ! Local variables

   logical ses 
   integer ip, iatm
   _REAL_ xi, yi, zi
   _REAL_ range1, rh
 
   epsx(1:xmymzm) = epsout; epsy(1:xmymzm) = epsout; epsz(1:xmymzm) = epsout

   if ( epsin /= epsout ) then

      atmsas(1:xmymzm) = 0

      ! use modified vdw surface (not ses) for cap water because there is only one
      ! sphere and also for the first FDPB because the grid spacing is too large

      if ( ifcap == 0 ) then
         if ( level == nfocus ) then
            ses = .true.
         else
            ses = .false.
         end if
      else
         ses = .false.
      end if
      rh = ONE/h

      ! solvent exluded surface

      if ( ses ) then

         ! part a, mark grid points within sa srf as 1, outside -1

         insas(1:xmymzm) = -1
         zv(1:xmymzm) = 9999.0d0
         do ip = 1, nsatm
            iatm = nzratm(ip); range1 = radip(iatm)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
            call exsasph( 1, insas, atmsas, zv(1) )
         end do
    
         ! save this for ion exclusion
    
         if (istrng /= ZERO) tv(1:xmymzm) = REAL(insas(1:xmymzm))
    
         ! part b, mark grid points within vdw srf as 2, outside 1 from part a
    
         zv(1:xmymzm) = 9999.0d0
         do ip = 1, nsatm
            iatm = nzratm(ip); range1 = radi(iatm)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
            call exvwsph( 2, insas, atmsas, zv(1) )
         end do
    
         ! part c, mark 1 grid points accessible to solvent probes as -2
    
         call contact( insas, atmsas )
    
         ! part d, mark grid points within solvent reentry probes as -1
    
         range1 = dprob*rh
         zv(1:xmymzm) = 9999.0d0
         do iatm = 1, narcdot
            xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
            call exresph(-1, insas, atmsas, zv(1) )
         end do
 
      ! modified vdw surface to approximate ses
 
      else
 
         ! mark volume within vdw srf as 2, outside -2, so there are only contact-type
         ! boundary grid points

         insas(1:xmymzm) = -2
         zv(1:xmymzm) = 9999.0d0
         do ip = 1, nsatm
            iatm = nzratm(ip);
            range1 = radip3(ip)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
            iatm = ip
            call exvwsph( 2, insas, atmsas, zv(1) )
         end do
 
      end if ! if ( ses ) then
 
      ! finally, save boundary edges for db energy and forces
 
      if ( level == nfocus ) call epsbnd( atmsas, insas )
 
      ! use the insas grid to setup epsx, epsy and epsz maps

      call epsmap( insas, atmsas, epsx, epsy, epsz )

   end if ! if ( epsin /= epsout ) then

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic sasurf
subroutine exsasph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   implicit none

   ! Passed variables

   integer  dielsph
   integer  insph(xm,ym,zm)
   integer  inatm(xm,ym,zm)
   _REAL_ dst(xm,ym,zm)
 
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d, d2, r

   r = rh*radi(iatm)

   lowk = ceiling(zi - range1); highk = floor(zi + range1)
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = ceiling(yi - range2); highj = floor(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = ceiling(xi - range3); highi = floor(xi + range3)
            do i = lowi, highi
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2; d = sqrt(d2)
               if ( d > r ) then
                  d = d - r
                  if ( insph(i,j,k) == dielsph ) then
                     if ( d < dst(i,j,k) ) then
                        inatm(i,j,k) = iatm; dst(i,j,k) = d
                     end if
                     cycle
                  end if
                  insph(i,j,k) = dielsph;
                  inatm(i,j,k) = iatm; dst(i,j,k) = d
               else
                  if ( insph(i,j,k) == dielsph ) cycle
                  insph(i,j,k) = dielsph; inatm(i,j,k) = iatm; dst(i,j,k) = ZERO
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
             
end subroutine exsasph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic vdw surf
subroutine exvwsph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   implicit none
    
   ! Passed variables
    
   integer  dielsph
   integer  insph(xm,ym,zm)
   integer  inatm(xm,ym,zm)
   _REAL_ dst(xm,ym,zm)
    
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d2
    
   lowk = ceiling(zi - range1); highk = floor(zi + range1)
   do k = lowk, highk
       
      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = ceiling(yi - range2); highj = floor(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = ceiling(xi - range3); highi = floor(xi + range3)
            do i = lowi, highi
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2
               if ( insph(i,j,k) == dielsph ) then
                  if ( d2 < dst(i,j,k) ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
                  cycle
               end if
               insph(i,j,k) = dielsph;
               inatm(i,j,k) = iatm; dst(i,j,k) = d2
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
    
end subroutine exvwsph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark contact grid points between vdw and sas surfaces
subroutine contact( insas,atmsas )
    
   implicit none
    
   integer insas(xm,ym,zm), atmsas(xm,ym,zm)
    
   integer i, j, k, buffer, iatm, ii, jj, iarc, inside
   _REAL_ xg(3), xi(3), xj(3), xij(3), dx(3), rxij, rd
   _REAL_ cosgij, cosgji, cosaij, cosaji

   _REAL_, parameter :: small = 0.01d0
    
   buffer = 1
   do k = 1+buffer, zm-buffer; do j = 1+buffer, ym-buffer; do i = 1+buffer, xm-buffer
       
      if ( insas(i,j,k) /= 1 ) cycle
       
      xg(1) = gox + i*h; xg(2) = goy + j*h; xg(3) = goz + k*h
       
      ! this is the atom that marked this grid within sasrf, so it will be the
      ! grid's conatct atom if it is marked so.
       
      iatm = atmsas(i,j,k)
      xi(1) = acrd(1,iatm); xi(2) = acrd(2,iatm); xi(3) = acrd(3,iatm)
       
      ! go through all arcs that this atom generates in circle() 
       
      inside = -2
      do ip = 1, marc(iatm)
         iarc = m2narc(ip,iatm)
          
         ! generated by outer loop, i.e. the atom is iatm in circle()
          
         if ( iarc >= fstarc(iatm) ) then
            jj = arcatm(1,iarc)
            xj(1) = acrd(1,jj)
            xj(2) = acrd(2,jj)
            xj(3) = acrd(3,jj)
            cosaij = savarc(1,iarc); cosaji = savarc(2,iarc)
          
         ! generated by inner loop, i.e. the atom is jatm in circle()
          
         else
            xj = xi
            ii = arcatm(2,iarc)
            xj(1) = acrd(1,ii)
            xj(2) = acrd(2,ii)
            xj(3) = acrd(3,ii)
            cosaji = savarc(1,iarc); cosaij = savarc(2,iarc)
         end if
         rxij = savarc(3,iarc)
         xij = rxij*(xj - xi)
          
         dx = xg - xi; rd = ONE/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))-small
          
         dx = xg - xj; rd = ONE/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))-small
          
         ! if gij < aij .and. gji < aji, this is a reentry grid
          
         if ( cosgij <= cosaij .or. cosgji <= cosaji ) cycle
         inside = 1
         exit
      end do
       
      insas(i,j,k) = inside
       
   end do; end do; end do
    
end subroutine contact
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within reentry surf
subroutine exresph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within a renentry sphere (dielectric constant dielsph)
   ! of index iatm as dielsph. Modified from UHBD (Comp. Phys. Comm. 91:57-95,
   ! 1995) routines excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! Passed variables

   integer  dielsph
   integer  insph(xm,ym,zm)
   integer  inatm(xm,ym,zm)
   _REAL_ dst(xm,ym,zm)

   ! Local variables

   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d2

   lowk = ceiling(zi - range1); highk = floor(zi + range1)
   do k = lowk, highk
       
      range2 = sqrt(range1**2-(zi-REAL(k))**2)
      lowj = ceiling(yi - range2); highj = floor(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-REAL(j))**2)
         if ( range3 > ZERO ) then
             
            lowi = ceiling(xi - range3); highi = floor(xi + range3)
            do i = lowi, highi
               if ( insph(i,j,k) == -2 ) cycle
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2
               if ( insph(i,j,k) == dielsph ) then
                  if ( d2 < dst(i,j,k) ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
                  cycle
               end if
               insph(i,j,k) = dielsph;
               inatm(i,j,k) = iatm; dst(i,j,k) = d2
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
    
end subroutine exresph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Save dielectric boundary grid points
subroutine epsbnd ( atmsas,insas )
    
   implicit none
    
   integer atmsas(xm,ym,zm), insas(xm,ym,zm)
    
   ! Local variables
    
   logical boundary
   integer buffer, i, j, k, clstmp
    
   nwarn = 0
   nbnd = 0
   buffer = 1
   do k = 1+buffer, zm-buffer; do j = 1+buffer, ym-buffer; do i = 1+buffer, xm-buffer
       
      ! set up condition for a boundary grid point
       
      boundary = .false.
      if ( (insas(i,j,k) ==  1 .or. insas(i,j,k) ==  2) .and.&
           (insas(i-1,j,k) == -1 .or. insas(i-1,j,k) == -2 .or. insas(i+1,j,k) == -1 .or.&
            insas(i+1,j,k) == -2 .or. insas(i,j-1,k) == -1 .or. insas(i,j-1,k) == -2 .or.&
            insas(i,j+1,k) == -1 .or. insas(i,j+1,k) == -2 .or. insas(i,j,k-1) == -1 .or.&
            insas(i,j,k-1) == -2 .or. insas(i,j,k+1) == -1 .or. insas(i,j,k+1) == -2) ) then 
            boundary = .true.
      else if ( (insas(i,j,k) == -1 .or. insas(i,j,k) == -2) .and.&
           (insas(i-1,j,k) ==  1 .or. insas(i-1,j,k) ==  2 .or. insas(i+1,j,k) ==  1 .or.&
            insas(i+1,j,k) ==  2 .or. insas(i,j-1,k) ==  1 .or. insas(i,j-1,k) ==  2 .or.&
            insas(i,j+1,k) ==  1 .or. insas(i,j+1,k) ==  2 .or. insas(i,j,k-1) ==  1 .or.&
            insas(i,j,k-1) ==  2 .or. insas(i,j,k+1) ==  1 .or. insas(i,j,k+1) ==  2) ) then
            boundary = .true.
      end if
      if ( .not. boundary ) cycle
 
      nbnd = nbnd + 1; iepsav(1,nbnd) = i; iepsav(2,nbnd) = j; iepsav(3,nbnd) = k
      if ( ifcap /= 0 ) then
         iepsav(4,nbnd) = -1
         cycle
      end if
 
      ! for a grid point in contact region +/- 2 or in a  solvent probe, simply use the atom/probe that
      ! marks it
 
      if ( abs(insas(i,j,k)) == 2 .or. insas(i,j,k) == -1 ) then
         clstmp = atmsas(i,j,k)
         if ( clstmp == 0 ) then
            if ( pbverbose ) write(6, '(a,4i4)') &
            'PB Warning in epsbnd(): No neighbor found for exposed boundary grid', i, j, k, insas(i,j,k)
            nwarn = nwarn + 1
         end if
 
      ! for a buried reentry grid point, find the atom that marked its neighoring exposed reentry
      ! grid points. Note that this may not be possible when grid spacing is large
 
      else if ( insas(i,j,k) == 1 ) then
         clstmp = fndcls( i, j, k, insas, atmsas )
         if ( clstmp == 0 ) then
            nwarn = nwarn + 1
         end if
      end if
 
      iepsav(4,nbnd) = clstmp
   end do; end do; end do
   if ( nwarn > 0 ) then
      if ( pbverbose ) write(6, '(a,i4)') &
      'PB Warning in epsbnd(): No neighbor found for boundary grids total:', nwarn
   end if
 
end subroutine epsbnd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find the closest reentry probe for a reentry boundary grid
function fndcls( i,j,k,insas,atmsas )
    
   implicit none
    
   ! Passed variables
    
   integer fndcls, i, j, k
   integer  insas(xm,ym,zm), atmsas(xm,ym,zm)
    
   ! Local variables
    
   integer iatm, l, lp, ip, jp, kp, iip(6), jjp(6), kkp(6), clsatm(6)
   _REAL_ xg, yg, zg
   _REAL_ dx, dy, dz, d, clsdst, clscrd(3,6)

   ! first stack these candidates into a 1-d list
    
   iip(1)=i-1; iip(2)=i+1; jjp(1:2)=j; kkp(1:2)=k
   iip(3:4)=i; jjp(3)=j-1; jjp(4)=j+1; kkp(3:4)=k
   iip(5:6)=i; jjp(5:6)=j; kkp(5)=k-1; kkp(6)=k+1
   lp = 0
   do l = 1, 6
      ip = iip(l); jp = jjp(l); kp = kkp(l)
      if ( atmsas(ip,jp,kp) == 0 .or. insas(ip,jp,kp) /= -1 ) cycle
      lp = lp + 1; iatm = atmsas(ip,jp,kp); clsatm(lp) = iatm
      clscrd(1,lp) = arccrd(1,iatm)
      clscrd(2,lp) = arccrd(2,iatm)
      clscrd(3,lp) = arccrd(3,iatm)
   end do
 
   ! now find the closest
 
   xg = gox + i*h; yg = goy + j*h; zg = goz + k*h
   clsdst = 999.d0
   fndcls = 0
   do ip = 1, lp
      dx = clscrd(1,ip) - xg; dy = clscrd(2,ip) - yg; dz = clscrd(3,ip) - zg
      d = abs(sqrt(dx**2 + dy**2 + dz**2) - dprob)
      if ( d >= clsdst ) cycle
      clsdst = d
      fndcls = clsatm(ip)
   end do
 
end function fndcls
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map insas into epsmap
subroutine epsmap( insas,atmsas,epsx,epsy,epsz )

   implicit none
   integer insas(xm,ym,zm), atmsas(xm,ym,zm)
   _REAL_ epsx(xm,ym,zm), epsy(xm,ym,zm), epsz(xm,ym,zm)

   integer i, j, k, a, b, c, d, a1, b1, c1, d1
   _REAL_ epsint

   epsint = TWO*epsin*epsout/(epsin+epsout)
 
   do k = 1, zm-1; do j = 1, ym-1; do i = 1, xm-1
      a = insas(i,j,k)
      b = insas(i+1,j,k)
      a1 = atmsas(i,j,k)
      b1 = atmsas(i+1,j,k)
      if ( sign(a,b) == a ) then
         if ( a > 0 ) then
            epsx(i,j,k) = epsin
         end if
      else
         if ( smoothopt == 1 ) call epsfracx(i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout)
         epsx(i,j,k) = epsint
      end if
      c = insas(i,j+1,k)
      c1 = atmsas(i,j+1,k)
      if ( sign(a,c) == a ) then
         if ( a > 0 ) then
            epsy(i,j,k) = epsin
         end if
      else
         if ( smoothopt == 1 ) call epsfracy(i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout)
         epsy(i,j,k) = epsint
      end if
      d = insas(i,j,k+1)
      d1 = atmsas(i,j,k+1)
      if ( sign(a,d) == a ) then
         if ( a > 0 ) then
            epsz(i,j,k) = epsin
         end if
      else
         if ( smoothopt == 1 ) call epsfracz(i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout)
         epsz(i,j,k) = epsint
      end if
   end do; end do; end do
 
!   do k = 1, zm
!      write(20, *) 'plane', k
!   do j = 1, ym
!      write(20, '(100f6.1)') epsx(1:xm,j,k)/eps0
!   end do
!   end do
!   do k = 1, zm
!      write(21, *) 'plane', k
!   do i = 1, xm
!      write(21, '(100f6.1)') epsy(i,1:ym,k)/eps0
!   end do
!   end do
!   do j = 1, ym
!      write(22, *) 'plane', j
!   do i = 1, xm
!      write(22, '(100f6.1)') epsz(i,j,1:zm)/eps0
!   end do
!   end do

end subroutine epsmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for x-edges
subroutine epsfracx( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! locate the atom that is crossing this x-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else if ( a == -1 ) then
      iatm = a1
   else if ( b == -1 ) then
      iatm = b1
   end if

   ! obtain the position and radius of the atom (probe) in grid unit

   if ( a == 2 .or. b == 2 ) then
      if ( ses ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh; iatm = nzratm(iatm)
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > ZERO ) then
         if ( b == 2 ) then
            aa = range3 - xi + REAL(i+1)
         else
            aa = range3 + xi - REAL(i)
         end if
         epsint = (depsout*depsin)/(depsin*(ONE-aa) + depsout*aa)
      else
         epsint = depsout
      end if
   else if ( a == -1 .or. b == -1 ) then
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > ZERO ) then
         if ( b == -1 ) then
            aa = range3 - xi + REAL(i+1)
         else
            aa = range3 + xi - REAL(i)
         end if
         epsint = (depsin*depsout)/(depsout*(ONE-aa) + depsin*aa)
      else
         epsint = depsin
      end if
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for y-edges
subroutine epsfracy( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! locate the atom that is crossing this y-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else if ( a == -1 ) then
      iatm = a1
   else if ( b == -1 ) then
      iatm = b1
   end if

   ! obtain the position and radius of the atom (probe) in grid unit

   if ( a == 2 .or. b == 2 ) then
      if ( ses ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh; iatm = nzratm(iatm)
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > ZERO ) then
         if ( b == 2 ) then
            aa = range3 - yi + REAL(j+1)
         else
            aa = range3 + yi - REAL(j)
         end if
         epsint = (depsout*depsin)/(depsin*(ONE-aa) + depsout*aa)
      else
         epsint = depsout
      end if
   else if ( a == -1 .or. b == -1 ) then
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > ZERO ) then
         if ( b == -1 ) then
            aa = range3 - yi + REAL(j+1)
         else
            aa = range3 + yi - REAL(j)
         end if
         epsint = (depsin*depsout)/(depsout*(ONE-aa) + depsin*aa)
      else
         epsint = depsin
      end if
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for z-edges
subroutine epsfracz( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! locate the atom that is crossing this z-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else if ( a == -1 ) then
      iatm = a1
   else if ( b == -1 ) then
      iatm = b1
   end if

   if ( a == 2 .or. b == 2 ) then
      if ( ses ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh; iatm = nzratm(iatm)
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
      if ( range3 > ZERO ) then
         if ( b == 2 ) then
            aa = range3 - zi + REAL(k+1)
         else
            aa = range3 + zi - REAL(k)
         end if
         epsint = (depsout*depsin)/(depsin*(ONE-aa) + depsout*aa)
      else
         epsint = depsout
      end if
   else if ( a == -1 .or. b == -1 ) then
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
      if ( range3 > ZERO ) then
         if ( b == -1 ) then
            aa = range3 - zi + REAL(k+1)
         else
            aa = range3 + zi - REAL(k)
         end if
         epsint = (depsin*depsout)/(depsout*(ONE-aa) + depsin*aa)
      else
         epsint = depsin
      end if
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracz

end subroutine pb_exmol
