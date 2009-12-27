
      ! inner loop over atoms...
      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and natom atomic charges

      eelrf = ZERO

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      eel = eel + eelrf

      ! now it is time to compute dbf
      ! first compute E_in and E_out ...
      !    compute E on the inner side of surface position crd(1:3)

      call gradu(xm,ym,zm,-ONE,8,4,fx0,fy0,fz0,dudxi0,dudyi0,dudzi0,phi,cphi,insas)

      dudxi0 = -dudxi0*FOURPI*eps0
      dudyi0 = -dudyi0*FOURPI*eps0
      dudzi0 = -dudzi0*FOURPI*eps0

 !    add the coulomb field to get the total E of inner side  

      dudxi0 = dudxi0 + coulomb(1)
      dudyi0 = dudyi0 + coulomb(2)
      dudzi0 = dudzi0 + coulomb(3)
      !    get the normal direction of position crd(1:3)
      !    then get the E of the outer side based on jump condition

      dn = sum(rn*rn)
      rn = rn / sqrt(dn)

      dudni = dudxi0*rn(1)+dudyi0*rn(2)+dudzi0*rn(3)
      dudno = dudni*epsin/epsout
      dudxo0 = dudxi0 - dudni*rn(1) + dudno*rn(1)
      dudyo0 = dudyi0 - dudni*rn(2) + dudno*rn(2)
      dudzo0 = dudzi0 - dudni*rn(3) + dudno*rn(3)
      EN = dudxi0*dudxo0+dudyi0*dudyo0+dudzi0*dudzo0

      if ( iatm > 0) then
         ds = h*h*dr/rsphere
         total_s = total_s + ds
      else
         ds = h*h*dr/rsphere
         total_s = total_s + ds
      end if

      f_bnd =  -INV_EIGHTPI*(epsout-epsin)*(dudni*dudni/epsout) ! for the normal field  approximation
!     f_bnd =  -INV_EIGHTPI*(epsout-epsin)*EN                 ! for the total field  
      dum(1) = f_bnd*rn(1)*ds
      dum(2) = f_bnd*rn(2)*ds
      dum(3) = f_bnd*rn(3)*ds

!     write(6,*) "ip, nbndx,iatm", ip, nbndx,iatm
      if ( iatm > 0  ) then
         f(1,iatm) = f(1,iatm) - dum(1)
         f(2,iatm) = f(2,iatm) - dum(2)
         f(3,iatm) = f(3,iatm) - dum(3)
      else
         iarc = dotarc(-iatm)
         matm = arcatm(1,iarc); natm = arcatm(2,iarc)
         d1=abs(sqrt((crd(1)-acrd(1,matm))**2+(crd(2)-acrd(2,matm))**2+(crd(3)-acrd(3,matm))**2)-radi(matm))
         d2=abs(sqrt((crd(1)-acrd(1,natm))**2+(crd(2)-acrd(2,natm))**2+(crd(3)-acrd(3,natm))**2)-radi(natm))
         f(1:3,matm) = f(1:3,matm)-d2/(d1+d2)*dum(1:3)
         f(1:3,natm) = f(1:3,natm)-d1/(d1+d2)*dum(1:3)
      end if
