
      dfx = factor1*dx
      dfy = factor1*dy
      dfz = factor1*dz
      
      ! collecting atomic forces

      if ( iatm > 0 ) then
         f(1,iatm) = f(1,iatm) + dfx
         f(2,iatm) = f(2,iatm) + dfy
         f(3,iatm) = f(3,iatm) + dfz
      else
         iarc = dotarc(-iatm)
         matm = arcatm(1,iarc); natm = arcatm(2,iarc)
         d1=abs(sqrt((crd(1)-acrd(1,matm))**2+(crd(2)-acrd(2,matm))**2+(crd(3)-acrd(3,matm))**2)-radi(matm))
         d2=abs(sqrt((crd(1)-acrd(1,natm))**2+(crd(2)-acrd(2,natm))**2+(crd(3)-acrd(3,natm))**2)-radi(natm))
         f(1,matm) = f(1,matm) - d2/(d1+d2)*dfx
         f(2,matm) = f(2,matm) - d2/(d1+d2)*dfy
         f(3,matm) = f(3,matm) - d2/(d1+d2)*dfz
         f(1,natm) = f(1,natm) - d1/(d1+d2)*dfx
         f(2,natm) = f(2,natm) - d1/(d1+d2)*dfy
         f(3,natm) = f(3,natm) - d1/(d1+d2)*dfz
      end if
