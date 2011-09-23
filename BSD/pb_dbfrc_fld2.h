
      ! inner loop over atoms...
      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and natom atomic charges

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      ! now it is time to compute dbf
      ! part a: compute E_in and E_out ...
      !    compute E on the inner side of surface position crd(1:3)

      call gradu(xm,ym,zm,-ONE,8,4,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,insas)

      dudxi0 = -dudxi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudyi0 = -dudyi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudzi0 = -dudzi0*FOURPI*eps0*AMBER_ELECTROSTATIC

      !    add the coulomb field to get the total E of inner side
      !    convert to displacement, D, for consistency with other methods

      dudxi0 = (dudxi0 + coulomb(1))*epsin/eps0
      dudyi0 = (dudyi0 + coulomb(2))*epsin/eps0
      dudzi0 = (dudzi0 + coulomb(3))*epsin/eps0

      !    get the normal direction of position crd(1:3)
      !    get normal field component, which is continuous
      !    get D of the outer side based on the jump conditions

      rn = rn/sqrt( sum(rn*rn) )

      dudni = dudxi0*rn(1)+dudyi0*rn(2)+dudzi0*rn(3)
      !dudno = dudni
      !dudxo0 = (dudxi0 - dudni*rn(1))/epsin*epsout + dudno*rn(1)
      !dudyo0 = (dudyi0 - dudni*rn(2))/epsin*epsout + dudno*rn(2)
      !dudzo0 = (dudzi0 - dudni*rn(3))/epsin*epsout + dudno*rn(3)

      ! part b: apply the normal field approximation
      !         or use the total field

      E2 = dudni*dudni
      !E2 = dudxi0*dudxo0 + dudyi0*dudyo0 + dudzi0*dudzo0

      ! part c: compute the surface force element

      fbnd = HALF*INV_FOURPI*(eps0*(epsin-epsout)/(epsin*epsout))*E2*ds
      dum(1) = fbnd*rn(1)
      dum(2) = fbnd*rn(2)
      dum(3) = fbnd*rn(3)

      if ( iatm > 0  ) then
         dbx(iatm) = dbx(iatm) + dum(1)
         dby(iatm) = dby(iatm) + dum(2)
         dbz(iatm) = dbz(iatm) + dum(3)
      else

         ! retrieve the two atoms that form the arc

         iarc = dotarc(-iatm)
         matm = arcatm(1,iarc); natm = arcatm(2,iarc)

         mvec(1:3) = x(1:3) - acrd(1:3,matm)
         nvec(1:3) = x(1:3) - acrd(1:3,natm)
         mvec = mvec*ONE/sqrt(mvec(1)**2 + mvec(2)**2 + mvec(3)**2)
         nvec = nvec*ONE/sqrt(nvec(1)**2 + nvec(2)**2 + nvec(3)**2)

         mxnv(1) = mvec(2)*nvec(3) - nvec(2)*mvec(3)
         mxnv(2) = nvec(1)*mvec(3) - mvec(1)*nvec(3)
         mxnv(3) = mvec(1)*nvec(2) - nvec(1)*mvec(2)
         mxnv = mxnv*ONE/sqrt(mxnv(1)**2 + mxnv(2)**2 + mxnv(3)**2)

         ! split dum() into tangent and normal directions wrt the plan of mvec/nvec

         dumnorm = dum(1)*mxnv(1) + dum(2)*mxnv(2) + dum(3)*mxnv(3)
         dum_norm = dumnorm*mxnv; dum_tang = dum - dum_norm

         ! further split dum_tangent into mvec and nvec directions

         mdotn = mvec(1)*nvec(1) + mvec(2)*nvec(2) + mvec(3)*nvec(3)
         rmdotn2 = ONE/(ONE - mdotn**2)
         fdotm = dum_tang(1)*mvec(1) + dum_tang(2)*mvec(2) + dum_tang(3)*mvec(3)
         fdotn = dum_tang(1)*nvec(1) + dum_tang(2)*nvec(2) + dum_tang(3)*nvec(3)
         if ( fdotm < ZERO .and. fdotn < ZERO) then
            mvec = -mvec; nvec = -nvec
         else if ( fdotm < ZERO ) then
            mvec = -mvec
            mdotn = -mdotn
         else if ( fdotn < ZERO ) then
            nvec = -nvec
            mdotn = -mdotn
         end if
         fdotm = abs(fdotm); fdotn = abs(fdotn)

         dfm = (fdotm - fdotn*mdotn)*rmdotn2
         dfn = (fdotn - fdotm*mdotn)*rmdotn2
         dbx(matm) = dbx(matm) + dfm*mvec(1) + HALF*dum_norm(1)
         dby(matm) = dby(matm) + dfm*mvec(2) + HALF*dum_norm(2)
         dbz(matm) = dbz(matm) + dfm*mvec(3) + HALF*dum_norm(3)
         dbx(natm) = dbx(natm) + dfn*nvec(1) + HALF*dum_norm(1)
         dby(natm) = dby(natm) + dfn*nvec(2) + HALF*dum_norm(2)
         dbz(natm) = dbz(natm) + dfn*nvec(3) + HALF*dum_norm(3)

         !d1=abs(sqrt((crd(1)-acrd(1,matm))**2+(crd(2)-acrd(2,matm))**2+(crd(3)-acrd(3,matm))**2)-radi(matm))
         !d2=abs(sqrt((crd(1)-acrd(1,natm))**2+(crd(2)-acrd(2,natm))**2+(crd(3)-acrd(3,natm))**2)-radi(natm))
         !dbx(matm) = dbx(matm) + d2/(d1+d2)*dum(1)
         !dby(matm) = dby(matm) + d2/(d1+d2)*dum(2)
         !dbz(matm) = dbz(matm) + d2/(d1+d2)*dum(3)
         !dbx(natm) = dbx(natm) + d1/(d1+d2)*dum(1)
         !dby(natm) = dby(natm) + d1/(d1+d2)*dum(2)
         !dbz(natm) = dbz(natm) + d1/(d1+d2)*dum(3)
      end if
