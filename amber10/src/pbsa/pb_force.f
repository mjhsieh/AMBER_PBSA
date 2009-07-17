! <compile=optimized>
#include "copyright.h"
#include "is_copyright.h"
#include "dprec.h"
#include "is_def.h"
#include "def_time.h"

module poisson_boltzmann

   implicit none

   ! PBMD parameters

   _REAL_, parameter :: pbkb   = 1.3807D-23 / 1.6606D-27 / (1.00D+12)**2 * (1.00D+10)**2
   _REAL_, parameter :: fioni  = 6.0220D+23 / 1.00D+30
   _REAL_, parameter :: fiono  = ONE / fioni
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   _REAL_, parameter :: frcfac = 0.01D0 / 4.1840D0

   ! PBMD FD control variables

   logical  :: outphi
   logical  :: srsas
   logical  :: scalerf
  
   integer  :: phiform 
   integer  :: dbfopt
   integer  :: solvopt
   integer  :: bcopt
   integer  :: xm
   integer  :: ym
   integer  :: zm
   integer  :: xmym
   integer  :: xmymzm
   integer  :: nbuffer
   integer  :: level
   integer  :: nfocus
   integer  :: fscale
   integer  :: maxitn
   integer  :: itn
   integer  :: savbcopt(MAXLEVEL)
   integer  :: savxm(MAXLEVEL)
   integer  :: savym(MAXLEVEL)
   integer  :: savzm(MAXLEVEL)
   integer  :: savxmym(MAXLEVEL)
   integer  :: savxmymzm(MAXLEVEL)

   _REAL_ :: h
   _REAL_ :: gox
   _REAL_ :: goy
   _REAL_ :: goz
   _REAL_ :: fmiccg
   _REAL_ :: accept
   _REAL_ :: norm
   _REAL_ :: inorm
   _REAL_ :: xmax
   _REAL_ :: xmin
   _REAL_ :: ymax
   _REAL_ :: ymin
   _REAL_ :: zmax
   _REAL_ :: zmin
   _REAL_ :: gxmax
   _REAL_ :: gxmin
   _REAL_ :: gymax
   _REAL_ :: gymin
   _REAL_ :: gzmax
   _REAL_ :: gzmin
   _REAL_ :: savxbox(MAXLEVEL)
   _REAL_ :: savybox(MAXLEVEL)
   _REAL_ :: savzbox(MAXLEVEL)
   _REAL_ :: cxbox(MAXLEVEL)
   _REAL_ :: cybox(MAXLEVEL)
   _REAL_ :: czbox(MAXLEVEL)
   _REAL_ :: savh(MAXLEVEL)
   _REAL_ :: savgox(MAXLEVEL)
   _REAL_ :: savgoy(MAXLEVEL)
   _REAL_ :: savgoz(MAXLEVEL)
   _REAL_ :: offx
   _REAL_ :: offy
   _REAL_ :: offz
   _REAL_ :: fillratio

   _REAL_ :: epsin
   _REAL_ :: epsout
   _REAL_ :: pbkappa
   _REAL_ :: istrng
   _REAL_ :: pbtemp
   _REAL_ :: totcrg
   _REAL_ :: totcrgp
   _REAL_ :: totcrgn

   _REAL_ :: pbgamma_int
   _REAL_ :: pbgamma_ext
 
   ! PBMD topology information

   integer               :: lastp

   integer , allocatable ::    icrd(:,:)
   integer , allocatable ::  grdcrg(:,:)
   _REAL_, allocatable :: qgrdcrg(:)
   _REAL_, allocatable ::    gcrd(:,:)
   _REAL_, allocatable ::    acrd(:,:)
   _REAL_, allocatable ::    acrg(:)
   _REAL_, allocatable ::    gcrg(:,:)
 
   ! PBMD nblist information

   integer               :: maxnbr
   integer               :: maxnba
   _REAL_              :: cutres, cutnb, cutfd, cutsa
 
   integer , allocatable ::   nshrt(:)
   integer , allocatable ::     nex(:)
   integer , allocatable ::     iex(:,:)
   integer , allocatable :: iprlong(:)
   integer , allocatable :: iprshrt(:)
   integer , allocatable ::  iar1pb(:,:)
   _REAL_, allocatable ::   cn1pb(:)
   _REAL_, allocatable ::   cn2pb(:)
   _REAL_, allocatable ::   cn3pb(:)

   ! PBMD cap water simulation information

   integer               :: mpopt
   integer               :: lmax
   integer               :: inatm
   integer               :: outatm
   _REAL_              :: sepbuf

   ! PBMD FD arrays for force/energy calculations
 
   integer               :: smoothopt
   integer               :: nbnd
   integer               :: nwarn
   integer , allocatable :: iepsav(:,:)
   integer , allocatable ::  insas(:)
   integer , allocatable :: atmsas(:)
   _REAL_, allocatable :: fedgex(:,:)
   _REAL_, allocatable :: fedgey(:,:)
   _REAL_, allocatable :: fedgez(:,:)
   integer , allocatable :: fatomx(:,:)
   integer , allocatable :: fatomy(:,:)
   integer , allocatable :: fatomz(:,:)
   _REAL_, allocatable ::    sbv(:)
   _REAL_, allocatable ::   epsx(:)
   _REAL_, allocatable ::   epsy(:)
   _REAL_, allocatable ::   epsz(:)
   _REAL_, allocatable ::    phi(:)
 
   ! PBMD FD solver working arrays
 
   _REAL_, allocatable ::  ad(:)
   _REAL_, allocatable ::  rd(:)
   _REAL_, allocatable :: am1(:)
   _REAL_, allocatable :: am2(:)
   _REAL_, allocatable :: am3(:)
   _REAL_, allocatable ::  bv(:)
   _REAL_, allocatable ::  zv(:)
   _REAL_, allocatable ::  pv(:)
   _REAL_, allocatable ::  tv(:)
   _REAL_, allocatable ::  xs(:)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of PBMD energy and forces
subroutine pb_force( natom,nres,ntypes,ipres,iac,ico,natex,cn1,cn2,cg,x,f,enb,eel,eelrf )
    
   use solvent_accessibility, only : dprob, radi, radip, radip2, radip3, nzratm, &
                                     sa_init, sa_driver, sa_free

   ! Common variables
    
#  include "pb_md.h"
#  include "md.h"
#  include "box.h"
    
   ! Passed variables
    
   integer natom, nres, ntypes, ipres(*), iac(*), ico(*), natex(*)
   _REAL_ cn1(*), cn2(*), cg(natom), x(3,natom), f(3,natom)
   _REAL_ enb, eel, eelrf
 
   ! Local variables

   integer iatm, proatm, atmfirst, atmlast
   integer atmind(natom), outflag(natom)
   _REAL_ acg(natom)
   _REAL_ pbcutcap, pbxcap, pbycap, pbzcap
   _REAL_ eelrffd, eelrfmp
   _REAL_ pbfrc(3,natom)
 
   enb = ZERO; eel = ZERO; eelrf = ZERO
   eelrffd = ZERO; eelrfmp = ZERO
   pbfrc = ZERO

   if ( ifcap /= 0) then
      pbcutcap = cutcap+TWO; pbxcap = xcap; pbycap = ycap; pbzcap = zcap 
      radi(-1) = pbcutcap; acrd(1,-1) = pbxcap; acrd(2,-1) = pbycap; acrd(3,-1) = pbzcap
      radip3(1) = radi(-1); nzratm(1) = -1
   end if
 
   ! split atoms into internal/external and update nblist
 
   call timer_start(TIME_PBLIST)
   if ( ntnbr == 1 .and. mpopt == 1 ) then
      call pb_atmpart(pbverbose,pbprint,natom,ibgwat,ienwat,inatm,outatm,ipres,outflag,&
              pbxcap,pbycap,pbzcap,pbcutcap,sepbuf,x)
   end if
   if ( ntnbr == 1 .and. cutres > ZERO ) then
      call pb_reslist(pbverbose,pbprint,maxnbr,natom,nres,ibgwat,ienwat,&
              ntypes,ipres,iac,ico,natex,nshrt,iar1pb,iprlong,cutres,x)
   end if
   if ( ntnba == 1 .and. max(cutnb,cutsa,cutfd) > ZERO ) then
      call pb_atmlist(pbverbose,pbprint,maxnba,natom,ntypes,iac,ico,natex,nshrt,nex,iex,&
              iar1pb,iprlong,iprshrt,cutnb,cutsa,cutfd,cn1,cn2,cn1pb,cn2pb,cn3pb,cg,x)
   end if
   call pb_atmconv(mpopt,natom,ibgwat,ienwat,outflag,atmind,ipres,x,cg,acg)
   if ( ntnbr == 1 ) ntnbr = 0
   if ( ntnba == 1 ) ntnba = 0
   call timer_stop(TIME_PBLIST)
 
   call timer_start(TIME_PBSETUP)
   if ( mpopt /=2 .and. pbgrid ) then
      call pb_setgrd(pbverbose,pbprint,pbinit,ifcap,natom,&
              pbxcap,pbycap,pbzcap,pbcutcap,acrd(1,1))
   end if
   call timer_stop(TIME_PBSETUP)

   ! compute grid-independent sas calculations for dielectric assignment
   ! when ifcap /= 0, no need to comptue sas

   call timer_start(TIME_PBSETUP)
   if ( srsas .and. ifcap == 0 ) then
      call sa_init(pbverbose,pbprint,natom,dprob,radi,radip,radip2)
      call sa_driver(pbverbose,pbprint,natom,dosas,ndosas,npbstep,nsaslag,&
              acrd(1,1),iar1pb(1,0),iprshrt,nex,iex)
   end if
   call timer_stop(TIME_PBSETUP)

   ! add FD reaction field energy/force

   call timer_start(TIME_PBFDFRC)
   if ( mpopt /= 2 .and. epsout /= epsin ) then
      call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,natom,pbfrc,eelrffd)
   end if
   call timer_stop(TIME_PBFDFRC)

   ! clean up for sas calculations

   call timer_start(TIME_PBSETUP)
   if ( srsas .and. ifcap == 0 ) then
      call sa_free( dosas,ndosas )
   end if
   call timer_stop(TIME_PBSETUP)

   ! add MP reaction field energy/forces when ifcap /= 0
 
   call timer_start(TIME_PBMP)
   if ( mpopt /= 0 .and. epsout /= epsin ) then
      if ( mpopt == 1 ) then      ! multipole expansion for boundary atoms
         atmfirst = inatm + 1
         atmlast  = natom
      else if ( mpopt == 2 ) then ! multipole expansion for all atoms
         atmfirst = 1
         atmlast  = natom
      end if
      call pb_mpfrc(natom,atmfirst,atmlast,lmax,pbcutcap,pbxcap,pbycap,pbzcap,&
              epsin,epsout,acrg,acrd(1,1),pbfrc,eelrfmp)
   end if
   call timer_stop(TIME_PBMP)
 
   ! add direct coulombic and nonbonded forces
    
   call timer_start(TIME_PBDIRECT)
   if ( ibgwat /= 0 ) then
      proatm = ipres(ibgwat) - 1
   else
      proatm = natom
   end if
   if ( cutnb == ZERO ) then
      call pb_directnocut(natom,proatm,ibgwat,ienwat,ntypes,iac,ico,nex,iex,cn1,cn2,acg,&
              acrd(1,1),pbfrc,eel,enb)
   else
      call pb_directwtcut(natom,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1),pbfrc,eel,enb)
   end if
   call timer_stop(TIME_PBDIRECT)

   ! returning:
   ! adding the nonbonded forces to the MD forces

   if ( dbfopt == 0 ) then 
      eel = eel + eelrffd + eelrfmp
      eelrf = ZERO
   else
      eelrf = eelrffd + eelrfmp
   end if

   do iatm = 1, natom
      f(1,iatm) = f(1,iatm) + pbfrc(1,iatm)
      f(2,iatm) = f(2,iatm) + pbfrc(2,iatm)
      f(3,iatm) = f(3,iatm) + pbfrc(3,iatm)
   end do

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ convert passed coordinates and charges to the internal format
subroutine pb_atmconv( mpopt,natom,ibgwat,ienwat,outflag,atmind,ipres,x,cg,acg )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! Passed variables
   
   integer mpopt, natom, ibgwat, ienwat
   integer outflag(natom), atmind(natom), ipres(*) 
   _REAL_ x(3,natom), cg(natom), acg(natom)
    
   ! Local variables
    
   integer i, j, ifirst, ilast, iatm, ires, num
    
   if ( mpopt == 1 ) then
       
      ! copy reordered coord/charge to private arrays for pb/mp
      ! protein atoms go into internal portion
       
      ifirst = 1; ilast = ipres(ibgwat-1)
      do iatm = ifirst, ilast
         acrd(1,iatm) = x(1,iatm); acrd(2,iatm) = x(2,iatm); acrd(3,iatm) = x(3,iatm)
         acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm); atmind(iatm) = iatm
      end do
       
      ! water atoms go into internal/external portion
       
      ifirst = ipres(ibgwat); ilast = natom
      i = ifirst; j =  inatm + 1
      do iatm = ifirst, ilast
         if ( outflag(iatm) == 0 ) then
            acrd(1,i   ) = x(1,iatm); acrd(2,i   ) = x(2,iatm); acrd(3,i   ) = x(3,iatm)
            acrg(i   ) = cg(iatm)/18.2223d0; acg(i   ) = cg(iatm); atmind(i   ) = iatm
            i = i + 1
         else
            acrd(1,j   ) = x(1,iatm); acrd(2,j   ) = x(2,iatm); acrd(3,j   ) = x(3,iatm);
            acrg(j   ) = cg(iatm)/18.2223d0; acg(j   ) = cg(iatm); atmind(j   ) = iatm
            j = j + 1
         end if
      end do
       
   else
       
      do iatm = 1, natom
         acrd(1,iatm) = x(1,iatm); acrd(2,iatm) = x(2,iatm); acrd(3,iatm) = x(3,iatm)
         acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm)
      end do
   end if
 
end subroutine pb_atmconv

end subroutine pb_force
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver for FDPB forces and energy
subroutine pb_fdfrc( pbverbose,pbprint,pbgrid,ifcap,natom,pbfrc,eelrf )

   ! passed variables
    
   logical pbverbose, pbprint, pbgrid  
   integer ifcap, natom
   _REAL_ eelrf, eelself, eelcoul, pbfrc(3,natom)!, fnet(3)
    
   ! local variables
    
   integer iatm, xsoffset 
   _REAL_ rh, fcrd(3,natom)
   _REAL_ aa, bb, cc, aa1, bb1, cc1
   _REAL_ bb1cc1, bb_cc1, bb1cc, bb_cc
    
   ! do fdpb calculations upto nfocus

   xsoffset = 1
   do level = 1, nfocus
       
      ! retrieving saved grid data into working variables

      bcopt = savbcopt(level)
      xm = savxm(level); ym = savym(level); zm = savzm(level)
      xmym = xm*ym; xmymzm = xmym*zm
      h = savh(level)
      gox = savgox(level); goy = savgoy(level); goz = savgoz(level)
       
      ! grid-unit version of the acrd array
       
      rh = ONE/h
      if ( ifcap /= 0 ) then 
         gcrd(1,-1) = (acrd(1,-1) - gox)*rh
         gcrd(2,-1) = (acrd(2,-1) - goy)*rh
         gcrd(3,-1) = (acrd(3,-1) - goz)*rh
      end if
      do iatm = 1, natom
         gcrd(1,iatm) = (acrd(1,iatm) - gox)*rh
         gcrd(2,iatm) = (acrd(2,iatm) - goy)*rh
         gcrd(3,iatm) = (acrd(3,iatm) - goz)*rh
      end do
      do iatm = 1, natom
         icrd(1,iatm) = floor(gcrd(1,iatm))
         icrd(2,iatm) = floor(gcrd(2,iatm))
         icrd(3,iatm) = floor(gcrd(3,iatm))
         fcrd(1,iatm) =  REAL(icrd(1,iatm))
         fcrd(2,iatm) =  REAL(icrd(2,iatm))
         fcrd(3,iatm) =  REAL(icrd(3,iatm))
      end do
      do iatm = 1, natom
         aa = gcrd(1,iatm) - fcrd(1,iatm)
         bb = gcrd(2,iatm) - fcrd(2,iatm)
         cc = gcrd(3,iatm) - fcrd(3,iatm)
         bb1 = ONE - bb; cc1 = ONE - cc
         aa  = acrg(iatm)*aa; aa1 = acrg(iatm) - aa
         bb1cc1 = bb1*cc1; bb_cc1 = bb *cc1
         bb1cc  = bb1*cc ; bb_cc  = bb *cc
         gcrg(1,iatm) = aa1*bb1cc1; gcrg(2,iatm) = aa *bb1cc1
         gcrg(3,iatm) = aa1*bb_cc1; gcrg(4,iatm) = aa *bb_cc1
         gcrg(5,iatm) = aa1*bb1cc ; gcrg(6,iatm) = aa *bb1cc
         gcrg(7,iatm) = aa1*bb_cc ; gcrg(8,iatm) = aa *bb_cc
      end do
       
      ! set up dielectric map
      ! when ifcap == 0, do dielectric map assignment everystep
      ! when ifcap /= 0, do dielectric map assignment once only when grid is set up

      call timer_start(TIME_PBEPS)
      if ( ifcap == 0 .or. pbgrid ) call pb_exmol( pbverbose,ifcap,natom )

      ! set up ion exclusion map

      if (istrng /= ZERO) then
         tv(1:xmymzm) = REAL(insas(1:xmymzm)); pv(1:xmymzm) = ZERO
         call pb_exion( pv(1), tv(1) )
      end if
      call timer_stop(TIME_PBEPS)
 
      ! now call fdpb driver

      call timer_start(TIME_PBSOLV)
      call pb_fddrv( 1,natom,phi,xs(xsoffset),sbv )
      call timer_stop(TIME_PBSOLV)

      ! if requested, print a summary when the grid is set up

      if ( pbverbose .and. pbprint ) then
         call pb_print( ifcap )
         write(6, *) '  Iterations required        :', itn
         write(6, *) '  Norm of the constant vector:', inorm
         write(6, *) '  Norm of the residual vector:', norm
         write(6, *) '  Convergence achieved       :', norm/inorm
         write(6, *)
      end if  !  pbverbose .and. pbprint
       
      ! if requested, output phi map
       
      if ( outphi .and. level == nfocus ) then
         if ( phiform == 0 ) then
            write(6,*) 'writing potential map in delphi format'
            open(64,file='pbsa.phi',form="unformatted")
            write(64) ' start of phimap    '
            write(64) ' potential', ' ------ AMBER PBSA ------ phimap in kT/e (0.593kcal/mol-e)  '
            write(64) real((frcfac/0.593d0)*phi(1:xmymzm))
            write(64) ' end of phimap  '
            write(64) real(1.0d0/h),real(cxbox(level)),real(cybox(level)),real(czbox(level)),xm
            close(64)
         elseif ( phiform == 1 ) then
            write(6,*) 'writing potential map in amber format'
            open(64,file='pbsa.phi',form="formatted")
            write(64,*) '# the following data is provided:'
            write(64,*) '# h, gox, goy, goz'
            write(64,*) '# xm, ym, zm'
            write(64,*) '# phi(1:xmymzm) in kcal/mol-e'
            write(64,*) '# mapping between (i,j,k) and phi index:'
            write(64,*) '# i + xm * ( j-1 + ym * ( k-1 ) )'
            write(64,*) '# grid coordinates: xg = gox + h*i; '
            write(64,*) '# yg = goy + h*j; zg = goz + h*k'
            write(64,*) h, gox, goy, goz
            write(64,*) xm, ym, zm
            write(64,*) frcfac*phi(1:xmymzm)
            close(64)
         endif
      end if  ! outphi .and. level == nfocus

      xsoffset = xsoffset + savxmymzm(level) + savxmym(level)

   end do  !  level = 1, nfocus

   ! compute fdfrc by the qE option

   if ( dbfopt == 0 ) then

      ! compute total qE energy and forces and delete self energy, note that self forces are zero

      call pb_qefrc( natom, eelrf, eelself, pbfrc, phi )

      ! delete fd grid energy and forces for all close pairs

      call pb_fdcoulomb( natom, eelcoul, pbfrc )

      eelrf = frcfac*(eelrf - eelself - eelcoul)
      !fnet = ZERO
      !do iatm = 1, natom
      !   fnet(1) = fnet(1) + pbfrc(1, iatm)
      !   fnet(2) = fnet(2) + pbfrc(2, iatm)
      !   fnet(3) = fnet(3) + pbfrc(3, iatm)
      !end do
      !fnet = fnet/dble(natom)
      !do iatm = 1, natom
      !   pbfrc(1,iatm) = frcfac*(pbfrc(1,iatm) - fnet(1))
      !   pbfrc(2,iatm) = frcfac*(pbfrc(2,iatm) - fnet(2))
      !   pbfrc(3,iatm) = frcfac*(pbfrc(3,iatm) - fnet(3))
      !end do
      pbfrc = frcfac*pbfrc

   ! compute fdfrc by db option

   else

      if ( epsin /= epsout ) then
         call pb_dbene( pbverbose,pbprint,natom,eelrf,pbfrc,insas,phi,sbv )
      end if

   end if
 
end subroutine pb_fdfrc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Exlusion of ions from protein interior, tv holds flags of sas surf.
subroutine pb_exion ( pv,tv )
 
   ! Passed variables
 
   _REAL_ pv(xm,ym,zm), tv(xm,ym,zm)
 
   ! Local variables
 
   integer i, j, k, buffer
   _REAL_ exclusion
 
   ! for InsightII display
   !open (unit=55, file='ions.dot')
   !write (55, '("DOTS")')
   buffer = 1
   do k = 1+buffer, zm-buffer; do j = 1+buffer, ym-buffer; do i = 1+buffer, xm-buffer
      exclusion = ZERO
      if ( tv(i-1,j,k) > ZERO .or. tv(i  ,j,k) > ZERO .or. tv(i+1,j,k) > ZERO .or. &
           tv(i,j-1,k) > ZERO .or. tv(i,j+1,k) > ZERO .or.&
           tv(i,j,k-1) > ZERO .or. tv(i,j,k+1) > ZERO ) then
         exclusion = SIX
         ! for InsightII display
         !g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
         !write (55,'(4(f8.3,2x))') g(1:3), 300.
      end if
      pv(i,j,k) = exclusion
   end do; end do; end do
   ! for InsightII display
   !close(55)

end subroutine pb_exion
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Finite-difference algorithm driver
subroutine pb_fddrv( atmfirst,atmlast,phi,xs,sbv )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Solving A * x = b, where A is dielectric/salt map, b is charge map, x is phi
   ! map.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! Passed variables
    
   integer atmfirst, atmlast
   _REAL_ phi(xmymzm), xs(xmymzm+xmym), sbv(xmymzm)
    
   ! Local variables
    
   integer l
   _REAL_ factor, factor1
    
   ! initialize A and b to zero, initialize the working arrays to zero,
   ! except RD, the reciprocal of D, to ONE. save final dielectric map
   ! for boundary force calculation since the AM() array will be overwritten.
    
   ad(     1-xmym  :xmymzm+xmym) = ZERO
   rd(     1-xmym  :0          ) = ONE
   rd(     xmymzm+1:xmymzm+xmym) = ONE
   am1(    1-xmym  :0          ) = epsout
   am1(    1       :xmymzm     ) = epsx(1:xmymzm)
   am1(    xmymzm+1:xmymzm+xmym) = ZERO
   am2(    1-xmym  :0          ) = epsout
   am2(    1       :xmymzm     ) = epsy(1:xmymzm)
   am2(    xmymzm+1:xmymzm+xmym) = ZERO
   am3(    1-xmym  :0          ) = epsout
   am3(    1       :xmymzm     ) = epsz(1:xmymzm)
   am3(    xmymzm+1:xmymzm+xmym) = ZERO
   bv(     1-xmym  :xmymzm+xmym) = ZERO
   zv(     1-xmym  :0          ) = ZERO
   zv(     xmymzm+1:xmymzm+xmym) = ZERO
   pv(     1-xmym  :0          ) = ZERO
   pv(     xmymzm+1:xmymzm+xmym) = ZERO
   tv(     1-xmym  :0          ) = ZERO
   tv(     xmymzm+1:xmymzm+xmym) = ZERO
    
   ! generate diagonal matrix elements, using the dielectric and ion exclusion maps

   if (istrng == ZERO) then
      do l = 1, xmymzm
         ad(l) = am1(l-1   ) + am1(l     ) + am2(l-xm  ) + am2(l     ) + am3(l-xmym) + am3(l     )
      end do
   else
      factor = epsout*(h*pbkappa)**2
      factor1 = factor/SIX
      do l = 1, xmymzm
         ad(l) = factor-pv(l)*factor1 + &
                 am1(l-1   ) + am1(l     ) + am2(l-xm  ) + am2(l     ) + am3(l-xmym) + am3(l     )
      end do
   end if
 
   ! set dielectric value extending outside the map to zero to avoid double counting
   ! since the boundary has been taken care by the charge map, bv().
   ! first the lower side, then the upper side
 
   am1(1-xmym:0) = ZERO
   am2(1-xmym:0) = ZERO
   am3(1-xmym:0) = ZERO
   call pb_setupper( am1(1), am2(1), am3(1) ) 

   ! place charges on the grid points and save them for induced charge calculations
 
   call pb_crggrd( atmfirst, atmlast, bv(1) )
   if ( dbfopt == 1 ) then
      factor = h/frcfac*(18.2223**2)/(epsin/eps0)
      sbv(1:xmymzm) = bv(1:xmymzm)*factor
   end if

   ! set the boundary condition, note except the first level, focusing is needed
 
   call pb_bndcnd( bv(1) )

   ! enter the core iteration routine

   if (solvopt == 1) then
      call pb_iccg( phi, xs )
   else
      write(6, *) 'PB bomb in pb_fddrv(): unknown solver'
      call mexit(6, 1)
   end if
 
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ set dielectric value extending outside the map to zero to avoid double counting
!+ This is the upper sides
subroutine pb_setupper( am1, am2, am3 )

   _REAL_ am1(xm,ym,zm), am2(xm,ym,zm), am3(xm,ym,zm)
   integer i, j, k

   do j = 1, ym; do k = 1, zm
      am1(xm,j,k) = ZERO
   end do; end do
   do i = 1, xm; do k = 1, zm
      am2(i,ym,k) = ZERO
   end do; end do
   do i = 1, xm; do j = 1, ym
      am3(i,j,zm) = ZERO
   end do; end do
 
end subroutine pb_setupper
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map atomic charges onto grid
subroutine pb_crggrd ( atmfirst, atmlast, bv )

   integer atmfirst, atmlast
   _REAL_ bv(xm,ym,zm)
 
   integer iatm, i, j, k
   _REAL_ rh 
   rh = ONE/h
   do iatm = atmfirst, atmlast
      i = icrd(1, iatm); j = icrd(2, iatm); k = icrd(3, iatm)
      bv(i  ,j  ,k  ) = bv(i  ,j  ,k  ) + gcrg(1, iatm)*rh
      bv(i+1,j  ,k  ) = bv(i+1,j  ,k  ) + gcrg(2, iatm)*rh
      bv(i  ,j+1,k  ) = bv(i  ,j+1,k  ) + gcrg(3, iatm)*rh
      bv(i+1,j+1,k  ) = bv(i+1,j+1,k  ) + gcrg(4, iatm)*rh
      bv(i  ,j  ,k+1) = bv(i  ,j  ,k+1) + gcrg(5, iatm)*rh
      bv(i+1,j  ,k+1) = bv(i+1,j  ,k+1) + gcrg(6, iatm)*rh
      bv(i  ,j+1,k+1) = bv(i  ,j+1,k+1) + gcrg(7, iatm)*rh
      bv(i+1,j+1,k+1) = bv(i+1,j+1,k+1) + gcrg(8, iatm)*rh
   end do

end subroutine pb_crggrd    
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign Debye-Huckel potential for the boundary charge grid
subroutine pb_bndcnd( bv )
   
   use constants, only: FOURPI 

   ! Common variables
    
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
     
   ! Passed variables
    
   _REAL_ bv(xm,ym,zm)
    
   ! Local variables
    
   integer i, j, k, iatm, ngrdcrg
   integer xmtmp, ymtmp, zmtmp, ix, iy, iz, itmp, jtmp, ktmp, idx, idy, idz
   _REAL_ htmp, goxtmp, goytmp, goztmp
   _REAL_ qtmp
   _REAL_ x, y, z
   _REAL_ xi, yi, zi, aa, bb, cc, aa1, bb1, cc1
   _REAL_ r, rinv
   _REAL_ factor
    
   factor = ONE/(FOURPI)
    
   ! bcopt = 1
   ! use zero potential for boundary grid, the boundary will be all solvent
    
   if ( bcopt == 1 ) then
      write(6, *) "PB bomb in pb_bndcnd(): zero BC not supported"
      call mexit(6, 1)
    
   ! bcopt = 2
   ! molecule dipolar debye-huckel contribution. the boundary will be all solvent.
    
   else if ( bcopt == 2 ) then
      write(6, *) "PB bomb in pb_bndcnd(): molecular dipolar BC not supported"
      call mexit(6, 1)
    
   ! bcopt = 3
   ! sum of residue dipolar debye-huckel contribution. the boundary will be all solvent.
    
   else if ( bcopt == 3 ) then
      write(6, *) "PB bomb in pb_bndcnd(): residue dipolar BC not supported"
      call mexit(6, 1)
    
   ! bcopt = 4
   ! sum of atom charge debye-huckel contribution. the boundary will be all solvent.
    
   else if ( bcopt == 4 ) then
       
      write(6, *) "PB bomb in pb_bndcnd(): atomic charge BC not supported"
      call mexit(6, 1)
       
   ! bcopt = 5
   ! sum of grid charge debye-huckel contribution. the boundary will be all solvent.
    
   else if ( bcopt == 5 ) then
       
      ! get grid-based charges
       
      ngrdcrg = 0
      do k = 1, zm; do j = 1, ym; do i = 1, xm
         if ( bv(i,j,k) == ZERO ) cycle
         ngrdcrg = ngrdcrg + 1
         grdcrg(1,ngrdcrg) = i
         grdcrg(2,ngrdcrg) = j
         grdcrg(3,ngrdcrg) = k
         qgrdcrg(ngrdcrg) = bv(i,j,k)*h
      end do; end do; end do
       
      factor = factor/h
      do iatm = 1, ngrdcrg
         itmp = grdcrg(1,iatm); jtmp = grdcrg(2,iatm); ktmp = grdcrg(3,iatm); qtmp = factor*qgrdcrg(iatm)
           
         ! k=0 and k=zm+1 faces
         do j = 1, ym; do i = 1, xm
            idx = abs(i-itmp); idy = abs(j-jtmp); idz = ktmp
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(REAL(idx**2 + idy**2 + idz**2))
               bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h*r))*qtmp/r
            end if
             
            idz = abs(zm+1-ktmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(REAL(idx**2 + idy**2 + idz**2))
               bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do
           
         ! j=0 and ym+1 faces
           
         do k = 1, zm; do i = 1, xm
            idx = abs(i-itmp); idy  = jtmp; idz  = abs(k-ktmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(REAL(idx**2 + idy**2 + idz**2))
               bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
              
            idy = abs(ym+1-jtmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(REAL(idx**2 + idy**2 + idz**2))
               bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do
      
         ! i=0 and i=xm+1 faces
      
         do k = 1, zm; do j = 1, ym
            idx = itmp; idy = abs(j-jtmp); idz = abs(k-ktmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(REAL(idx**2 + idy**2 + idz**2))
               bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
            
            idx = abs(xm+1-itmp)
            if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
               rinv = green(idx,idy,idz)
               bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(REAL(idx**2 + idy**2 + idz**2))
               bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do
      end do  !  iatm = 1, ngrdcrg
    
   ! bcopt = 0
   ! electrostatic focusing
    
   else if ( bcopt == 0 ) then
      xmtmp  = savxm(level-1) ; ymtmp  = savym(level-1) ; zmtmp  = savzm(level-1)
      htmp   = savh(level-1)
      goxtmp = savgox(level-1); goytmp = savgoy(level-1); goztmp = savgoz(level-1)
       
      ! k=0 and k=zm+1 faces
       
      do j = 1, ym; do i = 1, xm
          
         x  = gox + h*i        ; y  = goy + h*j        ; z  = goz
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa  = xi - REAL( ix ); bb  = yi - REAL( iy ); cc  = zi - REAL( iz )
         aa1 = ONE - aa       ; bb1 = ONE - bb       ; cc1 = ONE - cc
         bv(i,j,1 ) = bv(i,j,1 ) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         z  = goz + h*(zm+1)
         zi = (z - goztmp)/htmp
         iz = int( zi )
         cc  = zi - REAL( iz )
         cc1 = ONE - cc
         bv(i,j,zm) = bv(i,j,zm) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
   
      ! j=0 and j=ym+1 faces
   
      do k = 1, zm; do i = 1, xm
          
         x  = gox + h*i        ; y  = goy              ; z  = goz + h*k
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa = xi - REAL( ix ); bb = yi - REAL( iy ); cc = zi - REAL( iz )
         aa1 = ONE - aa      ; bb1 = ONE - bb      ; cc1 = ONE - cc
         bv(i,1 ,k) = bv(i,1 ,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         y  = goy + h*(ym+1)
         yi = (y - goytmp)/htmp
         iy = int( yi )
         bb  = yi - REAL( iy )
         bb1 = ONE - bb
         bv(i,ym,k) = bv(i,ym,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
        
      ! i=0 and i=xm+1 faces
        
      do k = 1, zm; do j = 1, ym
          
         x  = gox              ; y  = goy + h*j        ; z  = goz + h*k
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa  = xi - REAL( ix ); bb  = yi - REAL( iy ); cc  = zi - REAL( iz )
         aa1 = ONE - aa       ; bb1 = ONE - bb       ; cc1 = ONE - cc
         bv(1 ,j,k) = bv(1 ,j,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         x  = gox + h * (xm+1)
         xi = (x - goxtmp)/htmp
         ix = int( xi )
         aa  = xi - REAL( ix )
         aa1 = ONE - aa
         bv(xm,j,k) = bv(xm,j,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
       
   else
       
      ! unknown bcopt
       
      write(6, *) 'PB bomb in pb_bndcnd(): unknown BC option'
      call mexit(6, 1)
   end if  ! ( bcopt == 1 )
    
end subroutine pb_bndcnd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ phi interpretation, xm, ym, zm are for the previous phi map
_REAL_ function phintp( xmtmp,ymtmp,zmtmp,ix,iy,iz,aa,bb,cc,aa1,bb1,cc1 )
    
   ! Passed variables
    
   integer, intent(in) :: xmtmp, ymtmp, zmtmp, ix, iy, iz
   _REAL_, intent(in) :: aa, bb, cc, aa1, bb1, cc1
    
   ! Local Variables
    
   _REAL_ bb1cc1, bb_cc1, bb1cc, bb_cc
    
   ! determine the position of the point w.r.t. the map
    
   bb1cc1 = bb1*cc1; bb_cc1 = bb *cc1; bb1cc  = bb1*cc ; bb_cc  = bb *cc

   ! triliner interpolation
    
   phintp = aa1*bb1cc1*phi( ix   + xmtmp*( iy-1 + ymtmp*( iz-1 ) ) ) + &
            aa *bb1cc1*phi( ix+1 + xmtmp*( iy-1 + ymtmp*( iz-1 ) ) ) + &
            aa1*bb_cc1*phi( ix   + xmtmp*( iy   + ymtmp*( iz-1 ) ) ) + &
            aa *bb_cc1*phi( ix+1 + xmtmp*( iy   + ymtmp*( iz-1 ) ) ) + &
            aa1*bb1cc *phi( ix   + xmtmp*( iy-1 + ymtmp*( iz   ) ) ) + &
            aa *bb1cc *phi( ix+1 + xmtmp*( iy-1 + ymtmp*( iz   ) ) ) + &
            aa1*bb_cc *phi( ix   + xmtmp*( iy   + ymtmp*( iz   ) ) ) + &
            aa *bb_cc *phi( ix+1 + xmtmp*( iy   + ymtmp*( iz   ) ) )
                
end function phintp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ MICCG core iteration routine
subroutine pb_iccg ( phi, xs )
    
   ! Passed variables
    
   _REAL_ phi(xmymzm), xs(xmymzm+xmym)
    
   ! Local variables
    
   logical uconvg
   integer i, j
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2
    
   ! initialization
    
   do i = 1, xmymzm
      rd(i) = ONE/( ad(i) - &
         am1(i-1   )*(       am1(i-1   )+fmiccg*am2(i-1   )+fmiccg*am3(i-1   ))*rd(i-1   ) - &
         am2(i-xm  )*(fmiccg*am1(i-xm  )+       am2(i-xm  )+fmiccg*am3(i-xm  ))*rd(i-xm  ) - &
         am3(i-xmym)*(fmiccg*am1(i-xmym)+fmiccg*am2(i-xmym)+       am3(i-xmym))*rd(i-xmym) )
   end do
    
   do i = 1, xmymzm
      ad(i) = ad(i)*rd(i)
      rd(i) = sqrt(rd(i))
      bv(i) = bv(i)*rd(i)
      am1(i-1   ) = am1(i-1   )*rd(i)*rd(i-1   )
      am2(i-xm  ) = am2(i-xm  )*rd(i)*rd(i-xm  )
      am3(i-xmym) = am3(i-xmym)*rd(i)*rd(i-xmym)
   end do
    
   inorm = ZERO
   do i = 1, xmymzm
      ad(i) = ad(i) - TWO

      bv(i) = bv(i) + am1(i-1   )*bv(i-1   ) &
                    + am2(i-xm  )*bv(i-xm  ) &
                    + am3(i-xmym)*bv(i-xmym)
      inorm = inorm + abs(bv(i))
   end do
    
   do i = xmymzm, 1, -1
      tv(i) = xs(i) + am1(i     )*tv(i+1   ) &
                    + am2(i     )*tv(i+xm  ) &
                    + am3(i     )*tv(i+xmym)
   end do
   do i = 1, xmymzm
      zv(i) = xs(i) + ad (i     )*tv(i     ) &
                    + am1(i-1   )*zv(i-1   ) &
                    + am2(i-xm  )*zv(i-xm  ) &
                    + am3(i-xmym)*zv(i-xmym)
   end do
   bdotb1 = ZERO
   do i = xmymzm, 1, -1
      zv(i) = zv(i) + tv(i)
      bv(i) = bv(i) - zv(i)
      
      ! iteration 0.
      
      bdotb1 = bdotb1 + bv(i)*bv(i)
      pv(i)  = bv(i)
      
      ! first step of the matrix vector multiplication, see below
      
      tv(i) = pv(i) + am1(i     )*tv(i+1   ) &
                    + am2(i     )*tv(i+xm  ) &
                    + am3(i     )*tv(i+xmym)
   end do
    
   itn = 0
   uconvg = .true.
    
   ! the main loop of iccg solver
    
   do while ( uconvg )
       
      ! second and third steps of the matrix vector multiplication
       
      pdotz = ZERO
      do i = 1, xmymzm+xmym
         zv(i) = pv(i) + ad (i     )*tv(i     ) &
                       + am1(i-1   )*zv(i-1   ) &
                       + am2(i-xm  )*zv(i-xm  ) &
                       + am3(i-xmym)*zv(i-xmym)
         
         j = i - xmym
         zv(j) = zv(j) + tv(j)
         
         pdotz = pdotz + pv(j)*zv(j)
      end do
      alpha = bdotb1/pdotz
       
      norm = ZERO
      bdotb2 = ZERO
      itn = itn + 1
      do i = 1, xmymzm
         xs(i) = xs(i) + alpha*pv(i)
         bv(i) = bv(i) - alpha*zv(i)
         norm  = norm  + abs(bv(i))
         
         bdotb2= bdotb2+ bv(i)*bv(i)
      end do
       
      ! check convergence
       
      if ( itn >= maxitn .or. norm <= accept*inorm ) then
          
         uconvg = .false.
         if ( itn >= maxitn ) then
            write(6, *) 'PB warning in pb_miccg(): CG maxitn exceeded!'
         end if
          
      else
          
         beta = bdotb2/bdotb1
         bdotb1 = bdotb2
          
         ! first step of the matrix vector multiplication
          
         do i = xmymzm, 1, -1
            pv(i) = bv(i) + beta*pv(i)
             
            tv(i) = pv(i) + am1(i)*tv(i+1   ) &
                          + am2(i)*tv(i+xm  ) &
                          + am3(i)*tv(i+xmym)
         end do
      end if
   end do  !  while ( uconvg ), end of the main iccg loop
    
   ! back scaling of the solution
    
   do i = xmymzm, 1, -1
      tv(i)  = xs(i) + am1(i)*tv(i+1   ) &
                     + am2(i)*tv(i+xm  ) &
                     + am3(i)*tv(i+xmym)
       
      phi(i) = tv(i)*rd(i)
   end do
    
end subroutine pb_iccg

end subroutine pb_fddrv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ total finite difference es energy and forces
subroutine pb_qefrc( natom, grdnrg, grdself, pbfrc, phi )
   
   use constants, only: FOURPI 

   ! Common variables
   
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   
   ! Passed variables
   
   integer natom
   _REAL_ grdnrg, grdself
   _REAL_ pbfrc(3, natom)
   _REAL_ phi(xm,ym,zm)
   
   ! Local variables
   
   integer iatm
   integer i, j, k
   _REAL_ :: g000, g100, g110, g111
   _REAL_ :: gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   
   ! begin code
   
   g000 = green(0,0,0)/(FOURPI)
   g100 = green(1,0,0)/(FOURPI)
   g110 = green(1,1,0)/(FOURPI)
   g111 = green(1,1,1)/(FOURPI)
      
   grdnrg = ZERO
   grdself = ZERO
      
   ! split each atoms charge over the eight surrounding
   ! grid points according to the trilinear weighting
   ! function and add up each of the contributions.
      
   do iatm = 1, natom
      i = icrd(1,iatm); j = icrd(2,iatm); k = icrd(3,iatm)
      gci1 = gcrg(1,iatm); gci2 = gcrg(2,iatm)
      gci3 = gcrg(3,iatm); gci4 = gcrg(4,iatm)
      gci5 = gcrg(5,iatm); gci6 = gcrg(6,iatm)
      gci7 = gcrg(7,iatm); gci8 = gcrg(8,iatm)
      
      grdnrg = grdnrg + &
         gci1*phi(i  ,j  ,k  ) + gci2*phi(i+1,j  ,k  ) + &
         gci3*phi(i  ,j+1,k  ) + gci4*phi(i+1,j+1,k  ) + &
         gci5*phi(i  ,j  ,k+1) + gci6*phi(i+1,j  ,k+1) + &
         gci7*phi(i  ,j+1,k+1) + gci8*phi(i+1,j+1,k+1)
      
      grdself = grdself + &
         g000 * (gci1*gci1 + gci2*gci2 + gci3*gci3 + gci4*gci4 + &
                 gci5*gci5 + gci6*gci6 + gci7*gci7 + gci8*gci8 )*HALF + &
         g100 * (gci1*gci2 + gci1*gci3 + gci1*gci5 + gci2*gci4 + &
                 gci2*gci6 + gci4*gci3 + gci4*gci8 + gci3*gci7 + &
                 gci5*gci6 + gci5*gci7 + gci6*gci8 + gci8*gci7 ) + &
         g110 * (gci1*gci4 + gci1*gci6 + gci1*gci7 + gci2*gci3 + &
                 gci2*gci5 + gci2*gci8 + gci4*gci6 + gci4*gci7 + &
                 gci3*gci5 + gci3*gci8 + gci5*gci8 + gci6*gci7 ) + &
         g111 * (gci1*gci8 + gci2*gci7 + gci4*gci5 + gci6*gci3 )

      pbfrc(1,iatm) = &
         gci1*ex(i  ,j  ,k  ) + gci2*ex(i+1,j  ,k  ) + &
         gci3*ex(i  ,j+1,k  ) + gci4*ex(i+1,j+1,k  ) + &
         gci5*ex(i  ,j  ,k+1) + gci6*ex(i+1,j  ,k+1) + &
         gci7*ex(i  ,j+1,k+1) + gci8*ex(i+1,j+1,k+1)
      pbfrc(2,iatm) = &
         gci1*ey(i  ,j  ,k  ) + gci2*ey(i+1,j  ,k  ) + &
         gci3*ey(i  ,j+1,k  ) + gci4*ey(i+1,j+1,k  ) + &
         gci5*ey(i  ,j  ,k+1) + gci6*ey(i+1,j  ,k+1) + &
         gci7*ey(i  ,j+1,k+1) + gci8*ey(i+1,j+1,k+1) 
      pbfrc(3,iatm) = &
         gci1*ez(i  ,j  ,k  ) + gci2*ez(i+1,j  ,k  ) + &
         gci3*ez(i  ,j+1,k  ) + gci4*ez(i+1,j+1,k  ) + &
         gci5*ez(i  ,j  ,k+1) + gci6*ez(i+1,j  ,k+1) + &
         gci7*ez(i  ,j+1,k+1) + gci8*ez(i+1,j+1,k+1) 
   end do
   
   grdnrg  = HALF*grdnrg
   grdself = grdself/( h*epsin )

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ e field vector x
   _REAL_ function ex(i, j, k)

   integer, intent(in) :: i, j, k

   ex = ( phi(i-1,j  ,k  ) - phi(i+1,j  ,k  ) )/(2*h)

end function ex
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ e field vector y
_REAL_ function ey(i, j, k)

   integer, intent(in) :: i, j, k

   ey = ( phi(i  ,j-1,k  ) - phi(i  ,j+1,k  ) )/(2*h)

end function ey
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ e field vector z
_REAL_ function ez(i, j, k)

   integer, intent(in) :: i, j, k

   ez = ( phi(i  ,j  ,k-1) - phi(i  ,j  ,k+1) )/(2*h)

end function ez

end subroutine pb_qefrc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ FD coulombic energy and forces.
subroutine pb_fdcoulomb( natom, grdcoul, pbfrc )
   
   use constants, only: FOURPI 

   ! Common variables
   
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   
   ! Passed variables
   
   integer natom
   _REAL_ grdcoul
   _REAL_ pbfrc(3, natom)
   
   ! Local variables
   
   integer iatm, jatm, jp, ilast, jfirst, jlast
   integer dijx, dijx0, dijx1, dijx2, dijx3, dijy, dijy0, dijy1, dijy2, dijy3, dijz, dijz0, dijz1, dijz2, dijz3
   integer ix, iy, iz, jx, jy, jz
   _REAL_ gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ gcj1, gcj2, gcj3, gcj4, gcj5, gcj6, gcj7, gcj8
   _REAL_ gcij(27)
   _REAL_ factor, factor1
   _REAL_ ffx, ffy, ffz
   _REAL_ dumx, dumy, dumz
   _REAL_ frc(3, natom)
   
   ! begin code
   
   factor  = ONE/( FOURPI*epsin*h )
   factor1 = HALF/( FOURPI*epsin*h*h )
   
   grdcoul = ZERO
   frc = ZERO
   ilast = natom - 1
   do iatm = 1, ilast
      ix = icrd(1,iatm); iy = icrd(2,iatm); iz = icrd(3,iatm)
      
      gci1 = gcrg(1,iatm); gci2 = gcrg(2,iatm)
      gci3 = gcrg(3,iatm); gci4 = gcrg(4,iatm)
      gci5 = gcrg(5,iatm); gci6 = gcrg(6,iatm)
      gci7 = gcrg(7,iatm); gci8 = gcrg(8,iatm)
      
      jfirst = iar1pb(4, iatm-1) + 1
      jlast  = iar1pb(2, iatm)
      dumx = ZERO; dumy = ZERO; dumz = ZERO
      do jp = jfirst, jlast
         jatm = iprshrt(jp)
         jx = icrd(1,jatm); jy = icrd(2,jatm); jz = icrd(3,jatm)
         
         gcj1 = gcrg(1,jatm); gcj2 = gcrg(2,jatm)
         gcj3 = gcrg(3,jatm); gcj4 = gcrg(4,jatm)
         gcj5 = gcrg(5,jatm); gcj6 = gcrg(6,jatm)
         gcj7 = gcrg(7,jatm); gcj8 = gcrg(8,jatm)
         
         dijx  =       ix-jx; dijy  =       iy-jy; dijz  =       iz-jz
         dijx0 = abs(dijx-2); dijy0 = abs(dijy-2); dijz0 = abs(dijz-2)
         dijx1 = abs(dijx-1); dijy1 = abs(dijy-1); dijz1 = abs(dijz-1)
         dijx2 = abs(dijx+1); dijy2 = abs(dijy+1); dijz2 = abs(dijz+1)
         dijx3 = abs(dijx+2); dijy3 = abs(dijy+2); dijz3 = abs(dijz+2)
         dijx  = abs(dijx  ); dijy  = abs(dijy  ); dijz  = abs(dijz  )
         
         gcij( 1) = gci1*gcj1 + gci2*gcj2 + gci3*gcj3 + gci4*gcj4 + gci5*gcj5 + gci6*gcj6 + gci7*gcj7 + gci8*gcj8
         gcij( 2) = gci1*gcj2 + gci3*gcj4 + gci5*gcj6 + gci7*gcj8
         gcij( 3) = gci1*gcj3 + gci2*gcj4 + gci5*gcj7 + gci6*gcj8
         gcij( 4) = gci1*gcj4 + gci5*gcj8
         gcij( 5) = gci1*gcj5 + gci2*gcj6 + gci3*gcj7 + gci4*gcj8
         gcij( 6) = gci1*gcj6 + gci3*gcj8
         gcij( 7) = gci1*gcj7 + gci2*gcj8
         gcij( 8) = gci1*gcj8
         gcij( 9) = gci2*gcj1 + gci4*gcj3 + gci6*gcj5 + gci8*gcj7
         gcij(10) = gci2*gcj3 + gci6*gcj7
         gcij(11) = gci2*gcj5 + gci4*gcj7
         gcij(12) = gci2*gcj7
         gcij(13) = gci3*gcj1 + gci4*gcj2 + gci7*gcj5 + gci8*gcj6
         gcij(14) = gci3*gcj2 + gci7*gcj6
         gcij(15) = gci3*gcj5 + gci4*gcj6
         gcij(16) = gci3*gcj6
         gcij(17) = gci4*gcj1 + gci8*gcj5
         gcij(18) = gci4*gcj5
         gcij(19) = gci5*gcj1 + gci6*gcj2 + gci7*gcj3 + gci8*gcj4
         gcij(20) = gci5*gcj2 + gci7*gcj4
         gcij(21) = gci5*gcj3 + gci6*gcj4
         gcij(22) = gci5*gcj4
         gcij(23) = gci6*gcj1 + gci8*gcj3
         gcij(24) = gci6*gcj3
         gcij(25) = gci7*gcj1 + gci8*gcj2
         gcij(26) = gci7*gcj2
         gcij(27) = gci8*gcj1

         grdcoul = grdcoul + dble( &
              gci1*( gcj1*green(dijx ,dijy ,dijz ) + gcj2*green(dijx1,dijy ,dijz ) + &
                     gcj3*green(dijx ,dijy1,dijz ) + gcj4*green(dijx1,dijy1,dijz ) + &
                     gcj5*green(dijx ,dijy ,dijz1) + gcj6*green(dijx1,dijy ,dijz1) + &
                     gcj7*green(dijx ,dijy1,dijz1) + gcj8*green(dijx1,dijy1,dijz1) ) + &
              gci2*( gcj1*green(dijx2,dijy ,dijz ) + gcj2*green(dijx ,dijy ,dijz ) + &
                     gcj3*green(dijx2,dijy1,dijz ) + gcj4*green(dijx ,dijy1,dijz ) + &
                     gcj5*green(dijx2,dijy ,dijz1) + gcj6*green(dijx ,dijy ,dijz1) + &
                     gcj7*green(dijx2,dijy1,dijz1) + gcj8*green(dijx ,dijy1,dijz1) ) + &
              gci3*( gcj1*green(dijx ,dijy2,dijz ) + gcj2*green(dijx1,dijy2,dijz ) + &
                     gcj3*green(dijx ,dijy ,dijz ) + gcj4*green(dijx1,dijy ,dijz ) + &
                     gcj5*green(dijx ,dijy2,dijz1) + gcj6*green(dijx1,dijy2,dijz1) + &
                     gcj7*green(dijx ,dijy ,dijz1) + gcj8*green(dijx1,dijy ,dijz1) ) + &
              gci4*( gcj1*green(dijx2,dijy2,dijz ) + gcj2*green(dijx ,dijy2,dijz ) + &
                     gcj3*green(dijx2,dijy ,dijz ) + gcj4*green(dijx ,dijy ,dijz ) + &
                     gcj5*green(dijx2,dijy2,dijz1) + gcj6*green(dijx ,dijy2,dijz1) + &
                     gcj7*green(dijx2,dijy ,dijz1) + gcj8*green(dijx ,dijy ,dijz1) ) + &
              gci5*( gcj1*green(dijx ,dijy ,dijz2) + gcj2*green(dijx1,dijy ,dijz2) + &
                     gcj3*green(dijx ,dijy1,dijz2) + gcj4*green(dijx1,dijy1,dijz2) + &
                     gcj5*green(dijx ,dijy ,dijz ) + gcj6*green(dijx1,dijy ,dijz ) + &
                     gcj7*green(dijx ,dijy1,dijz ) + gcj8*green(dijx1,dijy1,dijz ) ) + &
              gci6*( gcj1*green(dijx2,dijy ,dijz2) + gcj2*green(dijx ,dijy ,dijz2) + &
                     gcj3*green(dijx2,dijy1,dijz2) + gcj4*green(dijx ,dijy1,dijz2) + &
                     gcj5*green(dijx2,dijy ,dijz ) + gcj6*green(dijx ,dijy ,dijz ) + &
                     gcj7*green(dijx2,dijy1,dijz ) + gcj8*green(dijx ,dijy1,dijz ) ) + &
              gci7*( gcj1*green(dijx ,dijy2,dijz2) + gcj2*green(dijx1,dijy2,dijz2) + &
                     gcj3*green(dijx ,dijy ,dijz2) + gcj4*green(dijx1,dijy ,dijz2) + &
                     gcj5*green(dijx ,dijy2,dijz ) + gcj6*green(dijx1,dijy2,dijz ) + &
                     gcj7*green(dijx ,dijy ,dijz ) + gcj8*green(dijx1,dijy ,dijz ) ) + &
              gci8*( gcj1*green(dijx2,dijy2,dijz2) + gcj2*green(dijx ,dijy2,dijz2) + &
                     gcj3*green(dijx2,dijy ,dijz2) + gcj4*green(dijx ,dijy ,dijz2) + &
                     gcj5*green(dijx2,dijy2,dijz ) + gcj6*green(dijx ,dijy2,dijz ) + &
                     gcj7*green(dijx2,dijy ,dijz ) + gcj8*green(dijx ,dijy ,dijz ) ) )
        
         ffx = dble( gcij( 1)*(green(dijx1,dijy ,dijz ) - green(dijx2,dijy ,dijz )) + &
                     gcij( 2)*(green(dijx0,dijy ,dijz ) - green(dijx ,dijy ,dijz )) + &
                     gcij( 3)*(green(dijx1,dijy1,dijz ) - green(dijx2,dijy1,dijz )) + &
                     gcij( 4)*(green(dijx0,dijy1,dijz ) - green(dijx ,dijy1,dijz )) + &
                     gcij( 5)*(green(dijx1,dijy ,dijz1) - green(dijx2,dijy ,dijz1)) + &
                     gcij( 6)*(green(dijx0,dijy ,dijz1) - green(dijx ,dijy ,dijz1)) + &
                     gcij( 7)*(green(dijx1,dijy1,dijz1) - green(dijx2,dijy1,dijz1)) + &
                     gcij( 8)*(green(dijx0,dijy1,dijz1) - green(dijx ,dijy1,dijz1)) + &
                     gcij( 9)*(green(dijx ,dijy ,dijz ) - green(dijx3,dijy ,dijz )) + &
                     gcij(10)*(green(dijx ,dijy1,dijz ) - green(dijx3,dijy1,dijz )) + &
                     gcij(11)*(green(dijx ,dijy ,dijz1) - green(dijx3,dijy ,dijz1)) + &
                     gcij(12)*(green(dijx ,dijy1,dijz1) - green(dijx3,dijy1,dijz1)) + &
                     gcij(13)*(green(dijx1,dijy2,dijz ) - green(dijx2,dijy2,dijz )) + &
                     gcij(14)*(green(dijx0,dijy2,dijz ) - green(dijx ,dijy2,dijz )) + &
                     gcij(15)*(green(dijx1,dijy2,dijz1) - green(dijx2,dijy2,dijz1)) + &
                     gcij(16)*(green(dijx0,dijy2,dijz1) - green(dijx ,dijy2,dijz1)) + &
                     gcij(17)*(green(dijx ,dijy2,dijz ) - green(dijx3,dijy2,dijz )) + &
                     gcij(18)*(green(dijx ,dijy2,dijz1) - green(dijx3,dijy2,dijz1)) + &
                     gcij(19)*(green(dijx1,dijy ,dijz2) - green(dijx2,dijy ,dijz2)) + &
                     gcij(20)*(green(dijx0,dijy ,dijz2) - green(dijx ,dijy ,dijz2)) + &
                     gcij(21)*(green(dijx1,dijy1,dijz2) - green(dijx2,dijy1,dijz2)) + &
                     gcij(22)*(green(dijx0,dijy1,dijz2) - green(dijx ,dijy1,dijz2)) + &
                     gcij(23)*(green(dijx ,dijy ,dijz2) - green(dijx3,dijy ,dijz2)) + &
                     gcij(24)*(green(dijx ,dijy1,dijz2) - green(dijx3,dijy1,dijz2)) + &
                     gcij(25)*(green(dijx1,dijy2,dijz2) - green(dijx2,dijy2,dijz2)) + &
                     gcij(26)*(green(dijx0,dijy2,dijz2) - green(dijx ,dijy2,dijz2)) + &
                     gcij(27)*(green(dijx ,dijy2,dijz2) - green(dijx3,dijy2,dijz2)) )
       
         ffy = dble( gcij( 1)*(green(dijx ,dijy1,dijz ) - green(dijx ,dijy2,dijz )) + &
                     gcij( 2)*(green(dijx1,dijy1,dijz ) - green(dijx1,dijy2,dijz )) + &
                     gcij( 3)*(green(dijx ,dijy0,dijz ) - green(dijx ,dijy ,dijz )) + &
                     gcij( 4)*(green(dijx1,dijy0,dijz ) - green(dijx1,dijy ,dijz )) + &
                     gcij( 5)*(green(dijx ,dijy1,dijz1) - green(dijx ,dijy2,dijz1)) + &
                     gcij( 6)*(green(dijx1,dijy1,dijz1) - green(dijx1,dijy2,dijz1)) + &
                     gcij( 7)*(green(dijx ,dijy0,dijz1) - green(dijx ,dijy ,dijz1)) + &
                     gcij( 8)*(green(dijx1,dijy0,dijz1) - green(dijx1,dijy ,dijz1)) + &
                     gcij( 9)*(green(dijx2,dijy1,dijz ) - green(dijx2,dijy2,dijz )) + &
                     gcij(10)*(green(dijx2,dijy0,dijz ) - green(dijx2,dijy ,dijz )) + &
                     gcij(11)*(green(dijx2,dijy1,dijz1) - green(dijx2,dijy2,dijz1)) + &
                     gcij(12)*(green(dijx2,dijy0,dijz1) - green(dijx2,dijy ,dijz1)) + &
                     gcij(13)*(green(dijx ,dijy ,dijz ) - green(dijx ,dijy3,dijz )) + &
                     gcij(14)*(green(dijx1,dijy ,dijz ) - green(dijx1,dijy3,dijz )) + &
                     gcij(15)*(green(dijx ,dijy ,dijz1) - green(dijx ,dijy3,dijz1)) + &
                     gcij(16)*(green(dijx1,dijy ,dijz1) - green(dijx1,dijy3,dijz1)) + &
                     gcij(17)*(green(dijx2,dijy ,dijz ) - green(dijx2,dijy3,dijz )) + &
                     gcij(18)*(green(dijx2,dijy ,dijz1) - green(dijx2,dijy3,dijz1)) + &
                     gcij(19)*(green(dijx ,dijy1,dijz2) - green(dijx ,dijy2,dijz2)) + &
                     gcij(20)*(green(dijx1,dijy1,dijz2) - green(dijx1,dijy2,dijz2)) + &
                     gcij(21)*(green(dijx ,dijy0,dijz2) - green(dijx ,dijy ,dijz2)) + &
                     gcij(22)*(green(dijx1,dijy0,dijz2) - green(dijx1,dijy ,dijz2)) + &
                     gcij(23)*(green(dijx2,dijy1,dijz2) - green(dijx2,dijy2,dijz2)) + &
                     gcij(24)*(green(dijx2,dijy0,dijz2) - green(dijx2,dijy ,dijz2)) + &
                     gcij(25)*(green(dijx ,dijy ,dijz2) - green(dijx ,dijy3,dijz2)) + &
                     gcij(26)*(green(dijx1,dijy ,dijz2) - green(dijx1,dijy3,dijz2)) + &
                     gcij(27)*(green(dijx2,dijy ,dijz2) - green(dijx2,dijy3,dijz2)) )
      
         ffz = dble( gcij( 1)*(green(dijx ,dijy ,dijz1) - green(dijx ,dijy ,dijz2)) + &
                     gcij( 2)*(green(dijx1,dijy ,dijz1) - green(dijx1,dijy ,dijz2)) + &
                     gcij( 3)*(green(dijx ,dijy1,dijz1) - green(dijx ,dijy1,dijz2)) + &
                     gcij( 4)*(green(dijx1,dijy1,dijz1) - green(dijx1,dijy1,dijz2)) + &
                     gcij( 5)*(green(dijx ,dijy ,dijz0) - green(dijx ,dijy ,dijz )) + &
                     gcij( 6)*(green(dijx1,dijy ,dijz0) - green(dijx1,dijy ,dijz )) + &
                     gcij( 7)*(green(dijx ,dijy1,dijz0) - green(dijx ,dijy1,dijz )) + &
                     gcij( 8)*(green(dijx1,dijy1,dijz0) - green(dijx1,dijy1,dijz )) + &
                     gcij( 9)*(green(dijx2,dijy ,dijz1) - green(dijx2,dijy ,dijz2)) + &
                     gcij(10)*(green(dijx2,dijy1,dijz1) - green(dijx2,dijy1,dijz2)) + &
                     gcij(11)*(green(dijx2,dijy ,dijz0) - green(dijx2,dijy ,dijz )) + &
                     gcij(12)*(green(dijx2,dijy1,dijz0) - green(dijx2,dijy1,dijz )) + &
                     gcij(13)*(green(dijx ,dijy2,dijz1) - green(dijx ,dijy2,dijz2)) + &
                     gcij(14)*(green(dijx1,dijy2,dijz1) - green(dijx1,dijy2,dijz2)) + &
                     gcij(15)*(green(dijx ,dijy2,dijz0) - green(dijx ,dijy2,dijz )) + &
                     gcij(16)*(green(dijx1,dijy2,dijz0) - green(dijx1,dijy2,dijz )) + &
                     gcij(17)*(green(dijx2,dijy2,dijz1) - green(dijx2,dijy2,dijz2)) + &
                     gcij(18)*(green(dijx2,dijy2,dijz0) - green(dijx2,dijy2,dijz )) + &
                     gcij(19)*(green(dijx ,dijy ,dijz ) - green(dijx ,dijy ,dijz3)) + &
                     gcij(20)*(green(dijx1,dijy ,dijz ) - green(dijx1,dijy ,dijz3)) + &
                     gcij(21)*(green(dijx ,dijy1,dijz ) - green(dijx ,dijy1,dijz3)) + &
                     gcij(22)*(green(dijx1,dijy1,dijz ) - green(dijx1,dijy1,dijz3)) + &
                     gcij(23)*(green(dijx2,dijy ,dijz ) - green(dijx2,dijy ,dijz3)) + &
                     gcij(24)*(green(dijx2,dijy1,dijz ) - green(dijx2,dijy1,dijz3)) + &
                     gcij(25)*(green(dijx ,dijy2,dijz ) - green(dijx ,dijy2,dijz3)) + &
                     gcij(26)*(green(dijx1,dijy2,dijz ) - green(dijx1,dijy2,dijz3)) + &
                     gcij(27)*(green(dijx2,dijy2,dijz ) - green(dijx2,dijy2,dijz3)) )
      
         dumx = dumx + ffx; dumy = dumy + ffy; dumz = dumz + ffz
         frc(1,jatm) = frc(1,jatm) - ffx
         frc(2,jatm) = frc(2,jatm) - ffy
         frc(3,jatm) = frc(3,jatm) - ffz
      end do  !  jp = jfirst, jlast
      
      frc(1,iatm) = frc(1,iatm) + dumx
      frc(2,iatm) = frc(2,iatm) + dumy
      frc(3,iatm) = frc(3,iatm) + dumz
   end do  !  iatm = 1, ilast
   
   grdcoul = factor*grdcoul
   pbfrc  = pbfrc - factor1*frc
   
end subroutine pb_fdcoulomb 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary energy and forces
subroutine pb_dbene( verbose,eneout,natom,eel,f,insas,phi,sbv )

   use constants, only: TWOPI, AMBER_ELECTROSTATIC2
   use solvent_accessibility, only : dprob, radi, arccrd, dotarc, arcatm
 
   ! Passed variables
 
   logical verbose, eneout
   integer natom
   integer insas(xm,ym,zm)
   _REAL_ phi(xm,ym,zm), sbv(xm,ym,zm)
   _REAL_ eel, f(3,natom)
 
   ! Local variables
 
   integer  i, j, k, iatm, jatm, matm, natm, iarc, ip
   _REAL_ srfcrg, factor, scalfact, r6
   _REAL_ g(3), x(3), dx(3), crd(3), dist, rdist, acg
   _REAL_ mvec(3), nvec(3), mdotn, mxnv(3), rmdotn2, fdotm, fdotn
   _REAL_ d2inv, dinv, de, dff, df(3), dfm, dfn, dum(3), dum_norm(3), dum_tang(3), dumnorm
   _REAL_ eelrf
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ fx(natom), fy(natom), fz(natom)

   _REAL_, parameter :: smallcrg = 0.5d0
 
   ! initialization
 
   r6 = SIXTH
   factor = THREE*h/(TWOPI)
   scalfact = ONE
 
   ax = acrd(1,1:natom)
   ay = acrd(2,1:natom)
   az = acrd(3,1:natom)

   srfcrg = ZERO; eel = ZERO
   fx = ZERO; fy = ZERO; fz = ZERO

   ! for InsightII display
   !open (unit=55, file='ms.dot')
   !write (55, '("DOTS")')
 
   do ip = 1, nbnd
      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k

      ! project the surface grid point on to the molecular surface, crd() is the new coord

      if (iatm == 0) then
         crd = g
      else
         if ( abs(insas(i,j,k)) == 2 ) then
            x(1:3) = acrd(1:3,iatm)
            dist = radi(iatm)
         else if ( abs(insas(i,j,k)) == 1 ) then
            x(1:3) = arccrd(1:3,iatm)
            dist = dprob
         end if
         dx = g - x
         rdist = dist*ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
         crd = x + dx*rdist
      end if
      ! for InsightII display
      !write (55,'(4(f8.3,2x))') crd(1:3), 300.
 
      ! compute induced charge on the molecular surface

      acg = factor*(phi(i,j,k)-&
      r6*( phi(i-1,j,k)+phi(i+1,j,k)+phi(i,j-1,k)+phi(i,j+1,k)+phi(i,j,k-1)+phi(i,j,k+1) ))
      acg = acg - sbv(i,j,k)
      srfcrg = srfcrg + acg*frcfac/(AMBER_ELECTROSTATIC2)
 
      ! compute reaction field energy and forces

      eelrf = ZERO
      dum = ZERO
      do jatm = 1, natom
         dx(1) = crd(1) - ax(jatm)
         dx(2) = crd(2) - ay(jatm)
         dx(3) = crd(3) - az(jatm)
         dinv = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2); d2inv = dinv**2

         de = acg*acrg(jatm)*dinv
         eelrf = eelrf + de

         dff = de*d2inv; df(1) = dx(1)*dff; df(2) = dx(2)*dff; df(3) = dx(3)*dff
         fx(jatm) = fx(jatm) + df(1)
         fy(jatm) = fy(jatm) + df(2)
         fz(jatm) = fz(jatm) + df(3)
         dum(1) = dum(1) - df(1)
         dum(2) = dum(2) - df(2)
         dum(3) = dum(3) - df(3)
      end do

      ! collecting energy

      eel = eel + eelrf
 
      ! collecting contact forces

      if ( abs(insas(i,j,k)) == 2 .and. iatm > 0 ) then
         fx(iatm) = fx(iatm) + dum(1)
         fy(iatm) = fy(iatm) + dum(2)
         fz(iatm) = fz(iatm) + dum(3)

      ! collecting reentry forces

      else if ( abs(insas(i,j,k)) == 1 .and. iatm > 0 ) then
         if ( iatm > 0 ) then
         iarc = dotarc(iatm)
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
         fx(matm) = fx(matm) + dfm*mvec(1) + HALF*dum_norm(1)
         fy(matm) = fy(matm) + dfm*mvec(2) + HALF*dum_norm(2)
         fz(matm) = fz(matm) + dfm*mvec(3) + HALF*dum_norm(3)
         fx(natm) = fx(natm) + dfn*nvec(1) + HALF*dum_norm(1)
         fy(natm) = fy(natm) + dfn*nvec(2) + HALF*dum_norm(2)
         fz(natm) = fz(natm) + dfn*nvec(3) + HALF*dum_norm(3)
         endif
      else
      end if

   end do

   ! for InsightII display
   !close(55)
   !stop

   if ( scalerf .and. abs(totcrg) > smallcrg ) then
      scalfact = abs( totcrg/srfcrg*(ONE/epsin - ONE/epsout)*eps0 )
      srfcrg = scalfact*srfcrg
      eel = scalfact*eel
      fx = scalfact*fx
      fy = scalfact*fy
      fz = scalfact*fz
   end if
    
   eel = HALF*frcfac*eel
   do iatm = 1, natom
      f(1,iatm) = f(1,iatm) - frcfac*fx(iatm)
      f(2,iatm) = f(2,iatm) - frcfac*fy(iatm)
      f(3,iatm) = f(3,iatm) - frcfac*fz(iatm)
   end do
    
   if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eel
   end if

end subroutine pb_dbene

end module poisson_boltzmann
