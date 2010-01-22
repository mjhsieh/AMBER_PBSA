! <compile=optimized>
#include "copyright.h"
#  define _REAL_ double precision
#include "pb_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for molecular ion map assignment
subroutine pb_ionmap( pbverbose,natom,iprob,&
              h,gox,goy,goz,xm,ym,zm,xmymzm,&
              gcrd,radi,&
              atmsas,insas,zv,iv)

   implicit none

   ! passed variables
 
   logical pbverbose
   integer natom
   _REAL_ iprob, h, gox, goy, goz
   integer xm, ym, zm, xmymzm
   _REAL_ gcrd(3,*), radi(*)
   integer atmsas(xmymzm), insas(xmymzm)
   _REAL_ zv(xmymzm), iv(xmymzm)

   ! local variables

   integer iatm
   _REAL_ rh, range0, range1, xi, yi, zi

   rh = 1.0d0/h

   ! resetting atmsas and insas

   insas(1:xmymzm) = -4; atmsas(1:xmymzm) = 0

   ! mark grid points within Stern layer as -3
    
   zv(1:xmymzm) = 9999.0d0
   do iatm = 1, natom
      range0 = radi(iatm)
      if ( range0 == 0.0d0 ) cycle
      range1 = (range0+iprob)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exstsph( -3, insas, atmsas, zv(1) )
   end do

   ! set up the ion exclusion map

   iv(1:xmymzm) = 1.0d0
   call ionmap( insas, iv )


contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic sasurf
subroutine exstsph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! Passed variables

   integer dielsph
   integer insph(xm,ym,zm)
   integer inatm(xm,ym,zm)
   _REAL_ dst(xm,ym,zm)
 
   ! Local variables
    
   integer i, j, k
   integer lowi, lowj, lowk
   integer highi, highj, highk
   _REAL_ range2, range3!, d, d2, r

   !r = range0

   lowk = ceiling(zi - range1); highk = floor(zi + range1)
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      lowj = ceiling(yi - range2); highj = floor(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then
             
            lowi = ceiling(xi - range3); highi = floor(xi + range3)
            do i = lowi, highi
               !d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2; d = sqrt(d2)
               !if ( d > r ) then
               !   d = d - r
               !   if ( insph(i,j,k) == dielsph ) then
               !      if ( d < dst(i,j,k) ) then
               !         inatm(i,j,k) = iatm; dst(i,j,k) = d
               !      end if
               !      cycle
               !   end if
               !   insph(i,j,k) = dielsph;
               !   inatm(i,j,k) = iatm; dst(i,j,k) = d
               !else
                  if ( insph(i,j,k) == dielsph ) cycle
                  insph(i,j,k) = dielsph; inatm(i,j,k) = iatm; dst(i,j,k) = 0.0d0
               !end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
             
end subroutine exstsph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Exlusion of ions from protein interior.
subroutine ionmap ( insas,stern )

   ! Passed variables

   integer insas(xm,ym,zm)
   _REAL_ stern(xm,ym,zm)

   ! Local variables

   integer i, j, k
   ! for InsightII display
   !_REAL_ g(3)

   ! for InsightII display
   !open (unit=55, file='ions.dot')
   !write (55, '("DOTS")')
   do k = 1, zm; do j = 1, ym; do i = 1, xm
      if ( insas(i,j,k) /= -4 ) then
         stern(i,j,k) = 0.0d0
         ! for InsightII display
         !g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
         !write (55,'(4(f8.3,2x))') g(1:3), 300.
      end if
   end do; end do; end do
   ! for InsightII display
   !close(55)

end subroutine ionmap


end subroutine pb_ionmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for molecular dielectric map assignment
subroutine pb_exmol_ses( pbverbose,ifcap,ipb,natom,&
              smoothopt,dprob,epsin,epsout,&
              h,gox,goy,goz,xm,ym,zm,xmymzm,level,nfocus,&
              nwarn,nsatm,narcdot,maxarc,nbnd,nbndx,nbndy,nbndz,&
              gcrd,acrd,radi,radip3,nzratm,&
              marc,m2narc,fstarc,arcatm,dotarc,arccrd,savarc,&
              atmsas,insas,lvlset,zv,epsx,epsy,epsz,&
              iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez,savbcopt )

   implicit none

   ! Passed variables
 
   logical pbverbose
   integer ifcap, ipb, natom, smoothopt
   _REAL_ dprob, epsin, epsout
   _REAL_ h, gox, goy, goz
   integer xm, ym, zm, xmymzm, level, nfocus
   integer nwarn, nsatm, narcdot, maxarc
   integer nbnd, nbndx, nbndy, nbndz
   _REAL_ gcrd(3,*), acrd(3,*), radi(*), radip3(*)
   integer nzratm(*), marc(*), m2narc(maxarc,*), fstarc(*), arcatm(2,*), dotarc(*)
   _REAL_ arccrd(3,*), savarc(3,*)
   integer atmsas(xmymzm), insas(xmymzm)
   _REAL_ lvlset(xmymzm), zv(xmymzm), epsx(xmymzm), epsy(xmymzm), epsz(xmymzm)
   integer iepsav(4,xmymzm), iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   _REAL_ fedgex(xmymzm), fedgey(xmymzm), fedgez(xmymzm)
   integer savbcopt(nfocus)
 
   ! Local variables

   integer ip, iatm, buf
   _REAL_ xi, yi, zi
   _REAL_ range1, rh
 
   epsx(1:xmymzm) = epsout; epsy(1:xmymzm) = epsout; epsz(1:xmymzm) = epsout

   ! in uniform dielectric systems, nothing to do here
   ! this is most likely for some reference state calculation or development work

   if ( epsin == epsout ) return

   ! in heterogeneous dielectrics systems ...

   rh = 1.0d0/h

   ! solvent exluded surface when at the fine level

   if ( level == nfocus ) then

      ! part a: reset atmsas and insas

      insas(1:xmymzm) = -4; atmsas(1:xmymzm) = 0

      ! part b: mark grid points just a bit larger than SAS by 2 grid points as -2

      ! WJ:
      ! 1. if use regularized SES, we first need to take care of the possibly very small
      ! solvent probe since we need level set function values on not just
      ! boundary grid points but its neighobrs as well. Adding 2 grids beyound
      ! dprob makes it safer for the later search of neighbors
      ! 2. this is also done for classical SES for comparison only ...

      zv(1:xmymzm) = 9999.0d0
      do iatm = 1, natom
         range1 = radi(iatm)
         if (range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh+2.0d0; xi = gcrd(1,iatm);yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( -2, insas, atmsas, zv(1) )
      end do

      ! part c, mark grid points within sa srf as 1, outside is -2 from b

      zv(1:xmymzm) = 9999.0d0
      do iatm = 1, natom
         range1 = radi(iatm)
         if ( range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( 1, insas, atmsas, zv(1) )
      end do
    
      ! part d, mark grid points within vdw srf as 2, outside is 1 from c
    
      zv(1:xmymzm) = 9999.0d0
      do iatm = 1, natom
         range1 = radi(iatm)
         if ( range1 == 0.0d0 ) cycle
         range1 = range1*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exvwsph( 2, insas, atmsas, zv(1) )
      end do
    
      ! part e, mark "1" grid points accessible to the solvent probe as -2
    
      call contact( insas, atmsas )
    
      ! part f, mark grid points within solvent reentry probes as -1
         
      buf = 2; range1 = dprob*rh + buf
      zv(1:xmymzm) = 9999.0d0
      do iatm = 1, narcdot
         xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
         call exresph( -1, insas, atmsas, zv(1) )
      end do
 

   ! modified vdw surface to approximate ses at the coarse level
 
   else
       
      ! part a, reset grid points everywhere as -2
 
      insas(1:xmymzm) = -2

      ! part b, mark volume within vdw srf as 2
      ! 1. since outside is -2, there are only contact fractional edges
      ! 2. note that zv is set to zero because we don't want to store any
      ! atom/distance info at this level, atmsas won't be used for this level
      ! 3. we are using the modified van der Waals radii that have been
      ! augamented by radinc in sa_driver()

      zv(1:xmymzm) = 0.0d0
      do ip = 1, nsatm
         iatm = nzratm(ip)
         range1 = radip3(ip)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exvwsph( 2, insas, atmsas, zv(1) )
      end do
 
   end if
 
   ! save boundary edges for eneopt=2 energy and forces
   ! Qin: always call epsbnd because bcopt == 6 requires this 
   ! WJ: epsbnd has been revised. should be free of warning in normal calls

   call epsbnd( atmsas, insas )

   ! WJ: if requested, set up level set function as signed distance to SES

   if ( ipb == 2 .and. level == nfocus ) then
      lvlset(1:xmymzm) = 9999.0d0
      call assignlvlset( atmsas, insas, lvlset )
   end if

   ! set up epsx, epsy and epsz maps
   ! boundary edges for eneopt=1 are saved inside.

   call epsmap( ipb, insas, atmsas, lvlset, epsx, epsy, epsz )


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

   lowk = max(1,ceiling(zi - range1)); highk = min(zm,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      lowj = max(1,ceiling(yi - range2)); highj = min(ym,floor(yi + range2))
      do j = lowj, highj

         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then

            lowi = max(1,ceiling(xi - range3)); highi = min(xm,floor(xi + range3))
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
                  insph(i,j,k) = dielsph; inatm(i,j,k) = iatm; dst(i,j,k) = 0.0d0
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
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
    
   lowk = max(1,ceiling(zi - range1)); highk = min(zm,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      lowj = max(1,ceiling(yi - range2)); highj = min(ym,floor(yi + range2))
      do j = lowj, highj

         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then

            lowi = max(1,ceiling(xi - range3)); highi = min(xm,floor(xi + range3))
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
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
   
 
end subroutine exvwsph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute density at grid points within mol sas surf
subroutine exdensph(ip,insph,packing,dens,atmctr )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Passed variables

   integer nn
   integer insph(xm,ym,zm)
   _REAL_ packing(xm,ym,zm), dens(xm,ym,zm)
   _REAL_ atmctr(xm,ym,zm)

   ! external function

   _REAL_ density

   ! Local variables
    
   integer  i, j, k, ii, jj, kk
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3,step
   _REAL_ d, r, d2, d2g, point, density_tmp
   integer  l , flag, ip
   
   step = 0.01d0
   d2g = 1.0d0/((2.0d0*dprob)*rh )

   lowk = int(zi - range1) + 1; highk = int(zi + range1)
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      lowj = int(yi - range2) + 1; highj = int(yi + range2)
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 > 0.0d0 ) then
             
            lowi = int(xi - range3) + 1; highi = int(xi + range3)
            do i = lowi, highi

               ! no need to spline if it is outside sas

               if ( insph(i,j,k) == 1) then
                  d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2; d = sqrt(d2)
                  d = d - r; point = d*d2g ! in the unit of (2*sprob/h)

                   ! get density by bicubic interpolation
                   
                   if( point > 1.0d0 ) cycle
                   if( point < -0.4 ) then
                      density_tmp = 3.352d0
                   else   
                      density_tmp = density( point, packing(i,j,k) )
                   end if
                   
                   dens(i,j,k) = dens(i,j,k) + density_tmp
                   !RLflag = i+xm*(j-1)+xmym*(k-1)
                   !RLif( density > density_value(1,flag) .and. density > density_value(2,flag)) then
                   !RL          density_list(2,flag)  = density_list(1,flag)   
                   !RL          density_list(1,flag)  = ip   
                   !RL          density_value(2,flag) = density_value(1,flag)   
                   !RL          density_value(1,flag) = density
                   !RLend if   
                   !RL          
                   !RLif( density < density_value(1,flag) .and. density > density_value(2,flag)) then
                   !RL          density_list(2,flag)  = ip   
                   !RL          density_value(2,flag) = density
                   !RLend if 
   !####################################################################################
   !               do ii =-40, 100
   !                 point = step*ii
   !                 if (abs(point) < 0.0000001) point = 0.0000001
   !                 call density_func(point,5.0d0,density, allc)
   !                 write(189,*)step*ii, density
   !               end do
   !               print *, 'OK-density' 
   !               stop 
                   
               end if ! ( insph(i,j,k) == 1 )

            end do  ! i = lowi, highi
             
         end if  ! ( range3 > ZERO )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
             
end subroutine exdensph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark contact grid points between vdw and sas surfaces
subroutine contact( insas,atmsas )
    
   integer insas(xm,ym,zm), atmsas(xm,ym,zm)
    
   integer i, j, k, buffer, iatm, ii, jj, iarc, inside
   _REAL_ xg(3), xi(3), xj(3), xij(3), dx(3), rxij, rd
   _REAL_ cosgij, cosgji, cosaij, cosaji

   _REAL_, parameter :: small = 0.01d0
    
   buffer = 1
   do k = 1+buffer, zm-buffer; do j = 1+buffer, ym-buffer; do i = 1+buffer, xm-buffer
       
      if ( insas(i,j,k) < 0 ) cycle
       
      xg(1) = gox + i*h; xg(2) = goy + j*h; xg(3) = goz + k*h
       
      ! this is the atom that marked this grid within sasrf, so it will be the
      ! grid's contact atom if it is marked so.
       
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
          
         dx = xg - xi; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))
          
         dx = xg - xj; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))
          
         ! if gij < aij .and. gji < aji, this is a reentry grid
          
         if ( cosgij <= cosaij .or. cosgji <= cosaji ) cycle
         if ( insas(i,j,k) == 1 ) then
            inside = 1
         else if ( insas(i,j,k) == 2 ) then
            inside = 1
         end if
         exit
      end do
       
      if ( inside == -2 .and. insas(i,j,k) == 2 ) inside = 2
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
       
      if ( k < 1 .or. k > zm ) cycle
      range2 = sqrt(range1**2-(zi-dble(k))**2)
      lowj = ceiling(yi - range2); highj = floor(yi + range2)
      do j = lowj, highj
          
         if ( j < 1 .or. j > ym ) cycle
         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then
             
            lowi = ceiling(xi - range3); highi = floor(xi + range3)
            do i = lowi, highi

               if ( i < 1 .or. i > xm ) cycle
               if ( insph(i,j,k) == -2 ) cycle
               if ( insph(i,j,k) == 2 ) cycle
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2
               if ( insph(i,j,k) == dielsph ) then
                  if ( d2 < dst(i,j,k) ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
                  cycle
               else if ( insph(i,j,k) == 1 ) then
                  if ( d2 < (range1 - buf)**2 ) then
                     insph(i,j,k) = dielsph;
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  else if ( d2 < dst(i,j,k) ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
               end if
!              insph(i,j,k) = dielsph;
!              inatm(i,j,k) = iatm; dst(i,j,k) = d2

            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
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
 
      clstmp = 0
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
         clstmp = atmsas(i,j,k)
         if ( clstmp == 0 ) then
            write(6,*) 'PB Info: Close contact cannot be found. fndcls() is called'
            clstmp = fndcls( i, j, k, insas, atmsas )
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
!+ Save dielectric boundary grid points
subroutine assignlvlset( atmsas,insas,u )

   implicit none

   ! passed variables

   integer atmsas(xm,ym,zm), insas(xm,ym,zm)
   _REAL_ u(xm,ym,zm)

   ! local variables

   integer i,j,k,l,buffer
   _REAL_ dist,d2,d

   rh = 1.0d0 / h
   dist = dprob*rh
   buffer = 1
   do l = 1, nbnd
      do k = iepsav(3,l) - buffer, iepsav(3,l) + buffer
         do j = iepsav(2,l) - buffer, iepsav(2,l) + buffer
            do i = iepsav(1,l) - buffer, iepsav(1,l) + buffer
               if ( atmsas(i,j,k) == 0 ) then
                  write(6,'(a,3i5)') 'PB Bomb in assignlvlset(): no atmsas', i,j,k
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l),iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l)+1,iepsav(2,l),iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l)-1,iepsav(2,l),iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l)+1,iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l)-1,iepsav(3,l))
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l),iepsav(3,l)+1)
                  write(6,*) iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l),iepsav(3,l)-1)
                  call mexit(6,1)
               end if
               if ( abs(insas(i,j,k)) == 2 ) then 
                  iatm = atmsas(i,j,k)
                  xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
                  d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2 ; d = sqrt(d2)
                  u(i,j,k) = dist - ( (radi(iatm)+dprob)*rh - d )
               else if ( abs(insas(i,j,k)) == 1 ) then
                  iatm = atmsas(i,j,k)
                  xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
                  d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2 ; d =sqrt(d2)
                  u(i,j,k) = dist - d
               else
                  write(6,'(a,4i5)') 'PB Bomb in assignlvlset(): illegal insas flag', i,j,k, insas(i,j,k)
                  call mexit(6,1)
               end if
            end do
         end do
      end do
   end do


end subroutine assignlvlset
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map insas into epsmap
subroutine epsmap( ipb,insas,atmsas,u,epsx,epsy,epsz )

   implicit none

   ! passed variables

   integer ipb
   integer insas(xm,ym,zm), atmsas(xm,ym,zm)
   _REAL_ u(xm,ym,zm)
   _REAL_ epsx(xm,ym,zm), epsy(xm,ym,zm), epsz(xm,ym,zm)

   ! local variables

   integer i, j, k, a, b, c, d, a1, b1, c1, d1
   integer x_flag, y_flag, z_flag
   _REAL_ epsint, epsint0

   ! set default value for simple harmonic average

   epsint0 = 2.0d0*epsin*epsout/(epsin+epsout)
   epsint = epsint0

   ! initialize boundary edge counters

   nbndx = 0; nbndy = 0; nbndz = 0
   x_flag = 0; y_flag = 0; z_flag = 0 ! local copies for the level set version

   ! the fraction caused by atom or probe

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
         if ( smoothopt > 0 .and. level == nfocus ) then
            if ( ipb == 1 .or. ipb == 2 .or. (ipb == 3 .and. (a == 2 .or. b == 2)) ) &
            call epsfracx  (i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout  )
            if ( (ipb == 2 .and. nbndx > x_flag) .or. (ipb == 3 .and. ( a == 1 .or. b == 1)) ) then
            call epsfracx_r(i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout,u)
            x_flag = x_flag + 1
            end if
            if ( smoothopt == 2 ) then
               if ( epsint > epsint0 ) then
                  epsint = epsout
               else
                  epsint = epsin
               end if
            end if
         end if
         epsx(i,j,k) = epsint
      end if
      c = insas(i,j+1,k)
      c1 = atmsas(i,j+1,k)
      if ( sign(a,c) == a ) then
         if ( a > 0 ) then
            epsy(i,j,k) = epsin
         end if
      else
         if ( smoothopt > 0 .and. level == nfocus ) then
            if ( ipb == 1 .or. ipb == 2 .or. (ipb == 3 .and. (a == 2 .or. c == 2)) ) &
            call epsfracy  (i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout  )
            if ( (ipb == 2 .and. nbndy > y_flag) .or. (ipb == 3 .and. ( a == 1 .or. c == 1)) ) then
            call epsfracy_r(i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout,u)
            y_flag = y_flag + 1
            end if
            if ( smoothopt == 2 ) then
               if ( epsint > epsint0 ) then
                  epsint = epsout
               else
                  epsint = epsin
               end if
            end if
         end if
         epsy(i,j,k) = epsint
      end if
      d = insas(i,j,k+1)
      d1 = atmsas(i,j,k+1)
      if ( sign(a,d) == a ) then
         if ( a > 0 ) then
            epsz(i,j,k) = epsin
         end if
      else
         if ( smoothopt > 0 .and. level == nfocus ) then
            if ( ipb == 1 .or. ipb == 2 .or. (ipb == 3 .and. (a == 2 .or. d == 2)) ) &
            call epsfracz  (i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout  )
            if ( (ipb == 2 .and. nbndz > z_flag) .or. (ipb == 3 .and. ( a == 1 .or. d == 1)) ) then
            call epsfracz_r(i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout,u)
            z_flag = z_flag + 1
            end if
            if ( smoothopt == 2 ) then
               if ( epsint > epsint0 ) then
                  epsint = epsout
               else
                  epsint = epsin
               end if
            end if
         end if
         epsz(i,j,k) = epsint
      end if

      ! checking the sanity of epsx, epsy, and epsz

      if ( epsx(i,j,k) < epsin .or. epsx(i,j,k) > epsout ) then
         write(6,'(a,3I10,a)') 'PB Bomb in epsmap(): epsx out of range', i,j,k
         call mexit(6,1)
      end if
      if ( epsy(i,j,k) < epsin .or. epsy(i,j,k) > epsout ) then
         write(6,'(a,3I10,a)') 'PB Bomb in epsmap(): epsy out of range', i,j,k
         call mexit(6,1)
      end if
      if ( epsz(i,j,k) < epsin .or. epsz(i,j,k) > epsout ) then
         write(6,'(a,3I10,a)') 'PB Bomb in epsmap(): epsz out of range', i,j,k
         call mexit(6,1)
      end if

   end do; end do; end do

   ! legacy checking output 
   !do k = 1, zm
   !   write(20, *) 'plane', k
   !do j = 1, ym
   !   write(20, '(100f6.1)') epsx(1:xm,j,k)/eps0
   !end do
   !end do
   !do k = 1, zm
   !   write(21, *) 'plane', k
   !do i = 1, xm
   !   write(21, '(100f6.1)') epsy(i,1:ym,k)/eps0
   !end do
   !end do
   !do j = 1, ym
   !   write(22, *) 'plane', j
   !do i = 1, xm
   !   write(22, '(100f6.1)') epsz(i,j,1:zm)/eps0
   !end do
   !end do

end subroutine epsmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for x-edges
subroutine epsfracx( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   integer flag, add_flag, flag_sub
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ front
   _REAL_ xg(3)

   ! locate the atom that is crossing this x-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else if ( a == -2 .and.b == 1 ) then
      iatm = a1
   else if ( b == -2 .and.a == 1 ) then
      iatm = b1
   else if ( a == 1  .and. b ==-1) then
      iatm = b1
   else if ( b == 1  .and. a ==-1) then
      iatm = a1
   else
      iatm = 0
   end if

   if ( iatm == 0 ) then
      write(6,*) 'PB Bomb in epsfrax(): no atom found to intersect an boundary edge'
      call mexit(6,1)
   end if

   ! always first obtain the position and radius of the atom (probe) in grid unit

   if ( a == 2 .or. b == 2 ) then
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh; iatm = nzratm(iatm)
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndx = nbndx + 1
         if ( b == 2 ) then
            aa = range3 - xi + dble(i+1)
            fedgex(nbndx) = 1.0d0 - aa
         else
            aa = range3 + xi - dble(i)
            fedgex(nbndx) = aa
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = iatm
      else
         epsint = depsout
      end if
   else if ( a == 1 .and. b == -2 ) then
      flag_sub = 1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(yi-j)**2
      if ( front >= 0.d0 ) then
         range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
         aa = xi+range3-dble(i)
         if ( .not. (aa < 0.d0 .or. aa > 1.d0 ) ) then 
            xg(1) = gox + h*i + h*aa
            xg(2) = goy + h*j
            xg(3) = goz + h*k
            call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
            if ( flag == 0 ) then ! the fraction does not lies in the reentry
               flag_sub = 0
               nbndx = nbndx + 1
               fedgex(nbndx) =  aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
         if ( range3 >= 0.0d0 ) then
            aa = xi-range3-dble(i)
            if ( .not. (aa < 0.d0 .or. aa > 1.d0) ) then
               flag_sub = 0
               nbndx = nbndx + 1
               fedgex(nbndx) =  aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = -iatm
            end if 
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( a == 1 .and. b == -1 ) then
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndx = nbndx + 1
         aa = xi-range3-dble(i)
         fedgex(nbndx) =  aa
         epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = -iatm
      else
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -2  ) then
      flag_sub = 1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(yi-j)**2
      if ( front >= 0.d0 ) then
         range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
         aa = dble(i+1)-xi+range3
         if ( .not. (aa < 0.d0 .or. aa > 1.d0 ) ) then 
            xg(1) = gox + h*i + h*(1-aa)
            xg(2) = goy + h*j
            xg(3) = goz + h*k
            call flag_value(xg,b1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndx = nbndx + 1
               fedgex(nbndx) = 1.0d0 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = b1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
         if ( range3 >= 0.0d0 ) then
            aa = dble(i+1)-xi-range3
            if ( .not. (aa < 0.d0 .or. aa > 1.d0) ) then
               flag_sub = 0 
               nbndx = nbndx + 1
               fedgex(nbndx) = 1.0d0 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -1  ) then
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh 
      yi = (arccrd(2,iatm) - goy)*rh 
      zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndx = nbndx + 1
         aa = dble(i+1)-xi-range3
         fedgex(nbndx) = 1.0d0 - aa
         epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = -iatm
      else
         epsint = depsin
      end if
   end if

   ! other situations will not be considered and use the default epsint value


end subroutine epsfracx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for y-edges
subroutine epsfracy( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm, flag
   integer add_flag, flag_sub
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ xg(3), front

   ! locate the atom that is crossing this y-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else if ( a == -2 .and. b == 1 ) then
      iatm = a1
   else if ( b == -2 .and. a == 1 ) then
      iatm = b1
   else if ( a == 1  .and. b ==-1 ) then
      iatm = b1
   else if ( b == 1  .and. a ==-1 ) then
      iatm = a1
   else
      iatm = 0
   end if

   if ( iatm == 0 ) then
      write(6,*) 'PB Bomb in epsfray(): no atom found to intersect an boundary edge'
      call mexit(6,1)
   end if

   ! always first obtain the position and radius of the atom (probe) in grid unit

   if ( a == 2 .or. b == 2 ) then
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh; iatm = nzratm(iatm)
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndy = nbndy + 1
         if ( b == 2 ) then
            aa = range3 - yi + dble(j+1)
            fedgey(nbndy) = 1.0d0 - aa
         else
            aa = range3 + yi - dble(j)
            fedgey(nbndy) = aa
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = iatm
      else
         epsint = depsout
      end if
   else if ( a == 1 .and. b == -2 ) then
      flag_sub = 1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(xi-i)**2
      if ( front >= 0.d0 ) then
         range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
         aa = yi+range3-dble(j)
         if ( .not. (aa < 0.d0 .or. aa > 1.d0 ) ) then 
            xg(1) = gox + h*i 
            xg(2) = goy + h*j + h*aa
            xg(3) = goz + h*k
            call flag_value(xg,a1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndy = nbndy + 1
               fedgey(nbndy) = aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = a1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
         if ( range3 >= 0.0d0 ) then
            aa = yi-range3-dble(j)
            if ( .not. (aa < 0.d0 .or. aa > 1.d0) ) then
               flag_sub = 0
               nbndy = nbndy + 1
               fedgey(nbndy) =  aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( a == 1 .and. b == -1 ) then
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh 
      yi = (arccrd(2,iatm) - goy)*rh 
      zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndy = nbndy + 1
         aa = yi-range3-dble(j)
         fedgey(nbndy) = aa
         epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = -iatm
      else
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -2 ) then
      flag_sub = 1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front =  range1**2-(zi-k)**2-(xi-i)**2
      if ( front >= 0.0d0 ) then
         range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
         aa = -yi+range3+dble(j+1)
         if ( .not. (aa < 0.d0 .or. aa > 1.d0 ) ) then 
            xg(1) = gox + h*i 
            xg(2) = goy + h*j + h*fedgey(nbndy)
            xg(3) = goz + h*k
            call flag_value(xg,b1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndy = nbndy + 1
               fedgey(nbndy) = 1.0d0 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = b1
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
         if ( range3 >= 0.0d0 ) then
            aa = dble(j+1) - yi - range3
            if ( .not. (aa < 0.d0 .or. aa > 1.d0) ) then
               flag_sub = 0 
               nbndy = nbndy + 1
               fedgey(nbndy) = 1.0d0 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -1 ) then
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh 
      yi = (arccrd(2,iatm) - goy)*rh 
      zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndy = nbndy + 1
         aa = -yi-range3+dble(j+1)
         fedgey(nbndy) = 1.0d0 - aa
         epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = -iatm
      else
         epsint = depsin
      end if
   end if

   ! other situations will not be considered and use the default epsint value


end subroutine epsfracy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for z-edges
subroutine epsfracz( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm, flag
   integer add_flag, flag_sub
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ xg(3), front

   ! locate the atom that is crossing this z-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else if ( a == -2 .and. b == 1 ) then
      iatm = a1
   else if ( b == -2 .and. a == 1 ) then
      iatm = b1
   else if ( a == 1  .and. b ==-1 ) then
      iatm = b1
   else if ( b == 1  .and. a ==-1 ) then
      iatm = a1
   else
      iatm = 0
   end if

   if ( iatm == 0 ) then
      write(6,*) 'PB Bomb in epsfraz(): no atom found to intersect an boundary edge'
      call mexit(6,1)
   end if

   ! always first obtain the position and radius of the atom (probe) in grid unit

   if ( a == 2 .or. b == 2 ) then
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh; iatm = nzratm(iatm)
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndz = nbndz + 1
         if ( b == 2 ) then
            aa = range3 - zi + dble(k+1)
            fedgez(nbndz) = 1.0d0 - aa
         else
            aa = range3 + zi - dble(k)
            fedgez(nbndz) = aa
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = iatm
      else
         epsint = depsout
      end if
   else if ( a == 1 .and. b == -2 ) then
      flag_sub = 1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(yi-j)**2-(xi-i)**2
      if ( front >= 0.d0 ) then
         range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
         aa = zi+range3-dble(k)
         if ( .not. (aa < 0.d0 .or. aa > 1.d0 ) ) then 
            xg(1) = gox + h*i 
            xg(2) = goy + h*j 
            xg(3) = goz + h*k + h*aa
            call flag_value(xg,a1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndz = nbndz + 1
               fedgez(nbndz) = aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = a1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
         if ( range3 >= 0.0d0 ) then
            aa = zi-range3-dble(k)
            if ( .not. (aa < 0.d0 .or. aa > 1.d0) ) then
               flag_sub = 0
               nbndz = nbndz + 1
               fedgez(nbndz) =  aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = -iatm
            end if 
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( a == 1 .and. b == -1 ) then
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh 
      yi = (arccrd(2,iatm) - goy)*rh 
      zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndz = nbndz + 1
         aa = zi-range3-dble(k)
         fedgez(nbndz) = aa
         epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = -iatm
      else
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -2 ) then
      flag_sub = 1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm) ; zi = gcrd(3,iatm)
      front = range1**2-(yi-j)**2-(xi-i)**2
      if ( front >= 0.0d0 ) then
         range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
         aa = -zi+range3+dble(k+1)
         if ( .not. (aa < 0.d0 .or. aa > 1.d0 ) ) then 
            xg(1) = gox + h*i 
            xg(2) = goy + h*j 
            xg(3) = goz + h*k + h*(1-aa)
            call flag_value(xg,b1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndz = nbndz + 1
               fedgez(nbndz) = 1.0d0 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = b1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
         if ( range3 >= 0.0d0 ) then
            aa = -zi+dble(k+1)-range3
            if ( .not. (aa < 0.d0 .or. aa > 1.d0) ) then
               flag_sub = 0
               nbndz = nbndz + 1 
               fedgez(nbndz) = 1.0d0 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -1 ) then
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh 
      yi = (arccrd(2,iatm) - goy)*rh 
      zi = (arccrd(3,iatm) - goz)*rh
      range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndz = nbndz + 1
         aa = -zi-range3+dble(k+1)
         fedgez(nbndz) = 1.0d0 - aa
         epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = -iatm
      else
         epsint = depsin
      end if
   end if

   ! other situations will not be considered and use the default epsint value


end subroutine epsfracz
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine flag_value(xg,iatm,flag)

   implicit none

   ! passed variables

   _REAL_ xg(3)
   integer iatm, flag

   ! local variables

   integer ii, jj, iarc
   _REAL_ xii(3), xjj(3), xij(3), dx(3), rxij, rd
   _REAL_ cosgij, cosgji, cosaij, cosaji, arcpos(3)
     
   arcpos(1:3) = arccrd(1:3,iatm)
   iarc = dotarc(iatm)

   ! generated by outer loop, i.e. the atom is iatm in circle()

   if ( iarc >= fstarc(iatm) ) then
      jj = arcatm(1,iarc)
      xjj(1) = acrd(1,jj)
      xjj(2) = acrd(2,jj)
      xjj(3) = acrd(3,jj)
      ii = arcatm(2,iarc)
      xii(1) = acrd(1,ii)
      xii(2) = acrd(2,ii)
      xii(3) = acrd(3,ii)
      cosaij = savarc(1,iarc); cosaji = savarc(2,iarc)

   ! generated by inner loop, i.e. the atom is jatm in circle()

   else
      ii = arcatm(2,iarc)
      xjj(1) = acrd(1,ii)
      xjj(2) = acrd(2,ii)
      xjj(3) = acrd(3,ii)
      jj = arcatm(1,iarc)
      xii(1) = acrd(1,jj)
      xii(2) = acrd(2,jj)
      xii(3) = acrd(3,jj)
      cosaji = savarc(1,iarc); cosaij = savarc(2,iarc)
   end if

   rxij = savarc(3,iarc)
   xij = rxij*(xjj - xii)

   dx = xg - xii; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
   cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

   dx = xg - xjj; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
   cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

   if ( cosgij > cosaij .and. cosgji > cosaji ) then
      flag = 1
   else
      flag = 0
   end if

end subroutine flag_value
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ assign fractional eps values for x-edges with the level set function
subroutine epsfracx_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout, u )

   implicit none

   ! passed variables

   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(xm,ym,zm)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0 

   if ( a > 0 ) then
      x1 = dble(i-1)
      x2 = dble(i  )
      x3 = dble(i+1)
      f1 = u(i-1,j,k)
      f2 = u(i  ,j,k)
      f3 = u(i+1,j,k)
      call root(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i+1,j,k) is inside, f2 < 0 and f1 > 0 

   if ( b > 0 ) then
      x1 = dble(i  )
      x2 = dble(i+1)
      x3 = dble(i+2)
      f1 = u(i  ,j,k)
      f2 = u(i+1,j,k)
      f3 = u(i+2,j,k)
      call root(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - real(i)
   else
      aa = real(i+1) - t
   end if

   fedgex(nbndx) = t - dble(i)
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)


end subroutine epsfracx_r
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ assign fractional eps values for y-edges with the level set function
subroutine epsfracy_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout, u )

   implicit none

   ! passed variables

   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(xm,ym,zm)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0 

   if ( a > 0 ) then
      x1 = dble(j-1)
      x2 = dble(j  )
      x3 = dble(j+1)
      f1 = u(i,j-1,k)
      f2 = u(i,j  ,k)
      f3 = u(i,j+1,k)
      call root(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i,j+1,k) is inside, f2 < 0 and f1 > 0

   if ( b > 0 ) then
      x1 = dble(j  )
      x2 = dble(j+1)
      x3 = dble(j+2)
      f1 = u(i,j  ,k)
      f2 = u(i,j+1,k)
      f3 = u(i,j+2,k)
      call root(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(j)
   else
      aa = dble(j+1) - t
   end if

   fedgey(nbndy) = t - dble(j)
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)

end subroutine epsfracy_r
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ assign fractional eps values for z-edges with the level set function
subroutine epsfracz_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout, u )

   implicit none

   ! passed variables

   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(xm,ym,zm)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0 

   if ( a > 0 ) then
      x1 = dble(k-1)
      x2 = dble(k)
      x3 = dble(k+1)
      f1 = u(i,j,k-1)
      f2 = u(i,j,k)
      f3 = u(i,j,k+1)
      call root(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i,j,k+1) is inside, f2 < 0 and f1 > 0

   if ( b > 0 ) then
      x1 = dble(k)
      x2 = dble(k+1)
      x3 = dble(k+2)
      f1 = u(i,j,k)
      f2 = u(i,j,k+1)
      f3 = u(i,j,k+2)
      call root(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(k)
   else
      aa = dble(k+1) - t
   end if

   fedgez(nbndz) =  t - dble(k)
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)


end subroutine epsfracz_r
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ root returns the approximated root between x0 and x1 if f0*f1 <=0 using
!+ quadratic interpolation
subroutine root(x0,x1,x2,f0,f1,f2,t0)

   implicit none

   ! passed variables

   _REAL_ x0,x1,x2,f0,f1,f2,t0

   ! local variables

   _REAL_ b,c,a0,b0,c0,t,r1,r2

   b = (f0-f1)/(x0-x1)
   c = f2 - f1 - b*(x2-x1)
   c = c/( (x2-x0)*(x2-x1))

   a0 = c
   b0 = b - c*(x0+x1)
   c0 = f1 -b*x1 + c*x0*x1

   if ( a0 == 0 ) then
      t0 = -c0/b0
      return
   end if

   t = b0*b0 - 4.0d0*a0*c0

   ! If t <=0, must be double root t is close to zero

   if ( t <= 0.0d0 ) then
      t0 = -b0/(2.0d0*a0)
      return
   end if

   t = sqrt(t)
   if ( b0 >= 0.0d0 ) then
      r1 = (-b0-t)/(2.0d0*a0)
   else
      r1 = (-b0+t)/(2.0d0*a0)
   end if

   r2 = -b0/a0-r1

   if ( x0 <= r1 + 1.0d-7 .and. r1 <= x1+1.0d-7 ) then
      t0 = r1
   else
      t0 = r2
   end if

   if ( x0 > t0 ) t0 = x0
   if ( x1 < t0 ) t0 = x1


end subroutine root
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute density function given distance and packing index
function density( distance,packing )

   implicit none

   ! common block variables
    
   _REAL_ allc(4,4,-2:4,0:14)
   common /density_coefficient/ allc

   ! passed variables

   _REAL_ density
   _REAL_ distance, packing

   ! local variables

   integer i, j, k, l
   _REAL_ u, v

   packing = 1.0d0 ! no need to use packing for now ...
 
   i = floor( distance/0.20d0 )
   j = floor( packing )

   if ( distance >= 0.0d0 ) then
      u = mod(distance/0.2d0, 1.0d0)
   else
      u = 1.0d0 - mod(abs(distance/0.2d0), 1.0d0)
   end if

   if ( packing >= 0.0d0 ) then
      v = mod(packing, 1.d0)
   else
      v = 1.0d0 - mod(abs(packing), 1.0d0)
   end if

   density = 0.0d0
   do l = 1, 4
      do k = 1, 4
         density = density + allc(k,l,i,j)*( u**(k-1) )*( v**(l-1) )
      end do
   end do


end function density


end subroutine pb_exmol_ses
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ set up coefficients for bicubic interpolation
subroutine density_init( )
    
   implicit none

   ! common block variables
    
   _REAL_ allc(4,4,-2:4,0:14)
   common /density_coefficient/ allc

   ! local variables

   integer i, j, k, l
   _REAL_ coef(4,4)
   _REAL_ f(4), f1(4), f12(4), f2(4)
   _REAL_ surf(-3:6,-1:16)
   _REAL_ surf1d(-2:5,-1:16)
   _REAL_ surf2d(-3:6,0:15)
   _REAL_ surf12d(-2:5,0:15)
   _REAL_ delta_dist, delta_pack

   !RLwrite(6,*) 'PB Info: Setting up spline coefficients for the density function'
   !RLopen(999,file='coef.dat')
   !RL
   !RL set up surf()
   !RL
   !RLdo j = -1, 16
   !RL   read(999,*) surf(-3:6,j)
   !RLend do
   !RL
   !RLclose(999)
 
   delta_dist = 0.20d0
   delta_pack = 1.0d0
 
   ! set up surf1d()
    
   do j = -1, 16
      do i = -2, 5
         surf1d(i,j) = (surf(i+1,j)-surf(i-1,j))/(delta_dist*2.0d0)
      end do
   end do
    
   ! set up surf2d()
    
   do j = 0, 15
      do i = -3, 6
         surf2d(i,j) = (surf(i,j+1)-surf(i,j-1))/(delta_pack*2.0d0)
      end do
   end do
    
   ! set up surf12d()
    
   do j = 0, 15
      do i = -2, 5
         surf12d(i,j) = (surf(i+1,j+1)-surf(i+1,j-1)-surf(i-1,j+1)+surf(i-1,j-1))/(delta_dist*delta_pack*4.0d0)
      end do
   end do
   
   ! now set up coefficients for each interpolation rectangle of x and y ...

   ! note the following order for surf(x, y): x changes first and y changes second, i.e. the fortran way
   ! both the orders for i,j and k,l are changed

   do j = 0, 14
      do i = -2, 4
          
         ! set up f
          
         f(1) = surf(i,j)
         f(2) = surf(i+1,j)
         f(3) = surf(i+1,j+1)
         f(4) = surf(i,j+1)
           
         ! set up f first derivitive of i
          
         f1(1) = surf1d(i,j)
         f1(2) = surf1d(i+1,j)
         f1(3) = surf1d(i+1,j+1)
         f1(4) = surf1d(i,j+1)
          
         ! set up f first derivitive of j
          
         f2(1) = surf2d(i,j)
         f2(2) = surf2d(i+1,j)
         f2(3) = surf2d(i+1,j+1)
         f2(4) = surf2d(i,j+1)
          
         ! set up f cross derivitive of i and j
         
         f12(1) = surf12d(i,j)
         f12(2) = surf12d(i+1,j)
         f12(3) = surf12d(i+1,j+1)
         f12(4) = surf12d(i,j+1)
         
         ! ready to call the set up routine ...
          
         call bicubic_coef(f,f1,f2,f12,delta_dist,delta_pack,coef)
          
         ! store the coefficients in the common block for later ...
          
         do l = 1,4
            do k = 1,4
               allc(l,k,i,j) = coef(l,k)
            end do
         end do
          
      end do
   end do


contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute coefficients for the bicubic interpolation
subroutine bicubic_coef(f,f1,f2,f12,d1,d2,coef)
   !
   ! the algorithm is documented in Numerical Recipes, 2nd edition, 1992.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! passed variables

   _REAL_ f(4), f1(4), f12(4), f2(4)
   _REAL_ d1, d2
   _REAL_ coef(4,4)

   ! local variables

   integer i, j, k
   _REAL_ d1d2, coef_tmp, ctmp(16), tmp(16)
   _REAL_ weight(16,16)

   ! weighting factors from Numerical Recipes

   !data wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,10*0, &
   !9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4, &
   !1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0, &
   !-6,4,2*0,3,-2,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2, &
   !10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4, &
   !-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0, &
   !2,-2,2*0,-1,1/
   data weight/ 1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4,&
                0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4,&
                0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4,&
                0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2,&
                0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2,&
                0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2,&
                0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2,&
                0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2,&
                0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1,&
                0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1,&
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1,&
                0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1 /

   ! setting up the working arrays

   d1d2 = d1*d2

   do i = 1, 4
      tmp(i) = f(i)
      tmp(i+4) = f1(i)*d1
      tmp(i+8) = f2(i)*d2
      tmp(i+12) = f12(i)*d1d2
   end do

   ! computing the coefficients

   do i = 1, 16
      ctmp(i) = 0.0d0
      do j = 1, 16
         ctmp(i) = ctmp(i) + weight(i,j)*tmp(j)
      end do
   end do

   ! saving it into the 2-d arragy for later

   k = 0
   do i = 1, 4
      do j = 1, 4
         k = k + 1
         coef(i,j) = ctmp(k)
      end do
   end do


end subroutine bicubic_coef


end subroutine density_init
