! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"
#include "timer.h"

module poisson_boltzmann

   implicit none

#  include "constants.h"

   ! PBMD parameters

   _REAL_, parameter :: pbkb   = 1.3807D-23 / 1.6606D-27 / (1.00D+12)**2 * (1.00D+10)**2
   _REAL_, parameter :: fioni  = 6.0220D+23 / 1.00D+30
   _REAL_, parameter :: fiono  = ONE / fioni
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   _REAL_, parameter :: frcfac = 0.01D0 / 4.1840D0

   ! PBMD FD control variables

   logical :: outphi
   logical :: srsas
   logical :: scalerf
  
   integer :: phiform 
   integer :: dbfopt
   integer :: eneopt
   integer :: frcopt
   integer :: frcdump
   integer :: solvopt
   integer :: npbopt
   integer :: bcopt
   integer :: smoothopt
   integer :: xm
   integer :: ym
   integer :: zm
   integer :: xmym
   integer :: xmymzm
   integer :: nbuffer
   integer :: nwarn
   integer :: level
   integer :: nfocus
   integer :: fscale
   integer :: maxitn
   integer :: itn
   integer :: savbcopt(MAXLEVEL)
   integer :: savxm(MAXLEVEL)
   integer :: savym(MAXLEVEL)
   integer :: savzm(MAXLEVEL)
   integer :: savxmym(MAXLEVEL)
   integer :: savxmymzm(MAXLEVEL)

   _REAL_ :: h
   _REAL_ :: gox
   _REAL_ :: goy
   _REAL_ :: goz
   _REAL_ :: fmiccg
   _REAL_ :: accept
   _REAL_ :: laccept
   _REAL_ :: wsor
   _REAL_ :: lwsor
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
   _REAL_ :: ivalence
   _REAL_ :: pbtemp
   _REAL_ :: totcrg
   _REAL_ :: totcrgp
   _REAL_ :: totcrgn

   _REAL_ :: pbgamma_int
   _REAL_ :: pbgamma_ext
 
   ! PBMD topology information

   integer :: lastp
   integer :: ngrdcrg

   integer, allocatable ::   icrd(:,:)
   integer, allocatable :: grdcrg(:,:)
   _REAL_, allocatable :: qgrdcrg(:)
   _REAL_, allocatable ::    gcrd(:,:)
   _REAL_, allocatable ::    acrd(:,:)
   _REAL_, allocatable ::    acrg(:)
   _REAL_, allocatable ::    gcrg(:,:)
 
   ! PBMD nblist information

   integer :: maxnbr
   integer :: maxnba
   _REAL_  :: cutres, cutnb, cutfd, cutsa
 
   integer, allocatable ::   nshrt(:)
   integer, allocatable ::     nex(:)
   integer, allocatable ::     iex(:,:)
   integer, allocatable :: iprlong(:)
   integer, allocatable :: iprshrt(:)
   integer, allocatable ::  iar1pb(:,:)
   _REAL_, allocatable ::    cn1pb(:)
   _REAL_, allocatable ::    cn2pb(:)
   _REAL_, allocatable ::    cn3pb(:)

   ! PBMD cap water simulation information

   integer :: mpopt
   integer :: lmax
   integer :: inatm
   integer :: outatm
   _REAL_  :: sepbuf

   ! physical variables for energy and force calculations

   integer :: nbnd
   integer :: nbndx
   integer :: nbndy
   integer :: nbndz

   ! physical variable maps for numerical solutions

   _REAL_, allocatable ::  phi(:)
   _REAL_, allocatable ::   bv(:)
   _REAL_, allocatable ::  sbv(:)
   _REAL_, allocatable :: epsx(:)
   _REAL_, allocatable :: epsy(:)
   _REAL_, allocatable :: epsz(:)
   _REAL_, allocatable ::   iv(:)

   ! geometry maps for dielectric interface
 
   integer, allocatable ::  insas(:)
   integer, allocatable :: atmsas(:)
   _REAL_, allocatable ::  lvlset(:)
   _REAL_, allocatable ::      zv(:)

   ! physical variable maps for force calculations

   _REAL_, allocatable ::     cphi(:)
   integer, allocatable ::  iepsav(:,:)
   integer, allocatable :: iepsavx(:,:)
   integer, allocatable :: iepsavz(:,:)
   integer, allocatable :: iepsavy(:,:)
   _REAL_, allocatable ::   fedgex(:)
   _REAL_, allocatable ::   fedgey(:)
   _REAL_, allocatable ::   fedgez(:)

   ! saved phi array for pbmd
 
   _REAL_, allocatable :: xs(:)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of PBMD energy and forces
subroutine pb_force( natom,nres,ntypes,ipres,iac,ico,natex,cn1,cn2,cg,x,f,enb,eel,eelrf )
    
   use solvent_accessibility, only : dprob, iprob, radi, radip, radip2, radip3, nzratm, &
       narcdot,maxarc,marc,m2narc,fstarc,arcatm,arccrd,savarc, &
       sa_init, sa_driver, sa_free
   use timer_module

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
   _REAL_ ionene
!  _REAL_ charge_ratio
 
!  open(866,file='coulomb_force.dat')

   enb = ZERO; eel = ZERO; eelrf = ZERO
   eelrffd = ZERO; eelrfmp = ZERO
   pbfrc = ZERO; ionene = ZERO
   atmind = 0

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
      call sa_driver(pbverbose,pbprint,ipb,inp,natom,dosas,ndosas,npbstep,nsaslag,&
              acrd(1,1),iar1pb(1,0),iprshrt,nex,iex)
   end if
   call timer_stop(TIME_PBSETUP)

   ! add FD reaction field energy/force

   call timer_start(TIME_PBFDFRC)
   if ( mpopt /= 2 .and. epsout /= epsin ) then
      call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,natom,pbfrc,eelrffd,ionene,npbstep,npbgrid,nstlim)
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

   pbfrc = ZERO ! resetting to zero for printing only
   if ( cutnb == ZERO ) then
      call pb_directnocut(natom,proatm,ibgwat,ienwat,ntypes,iac,ico,nex,iex,cn1,cn2,acg,&
              acrd(1,1),pbfrc,eel,enb)
   else
      call pb_directwtcut(natom,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1),pbfrc,eel,enb)
   end if
   call timer_stop(TIME_PBDIRECT)

   ! returning:
   ! adding the nonbonded forces to the MD forces

   if ( eneopt == 1 .and. savbcopt(1) /=6 ) then 
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
subroutine pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,natom,pbfrc,eelrf,ionene,npbstep,npbgrid,nstlim )

   use timer_module
   use solvent_accessibility, only : dprob,iprob,radi,radip3,nzratm,nsatm, &
       narcdot,maxarc,marc,m2narc,fstarc,arcatm,arccrd,savarc,dotarc

   implicit none

   ! passed variables
    
   logical pbverbose, pbprint, pbgrid  
   integer ifcap, ipb, natom, npbstep, npbgrid, nstlim
   _REAL_ ionene, eelrf, pbfrc(3,natom)!, fnet(3)
    
   ! local variables
    
   integer iatm, xsoffset 
   _REAL_ eelself, eelcoul
   _REAL_ rh, fcrd(3,natom)
   _REAL_ aa, bb, cc, aa1, bb1, cc1
   _REAL_ bb1cc1, bb_cc1, bb1cc, bb_cc
    
   ! do fdpb calculations upto nfocus

   eelself = ZERO; eelcoul = ZERO
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
      if ( ifcap == 0 .or. pbgrid ) then

         ! part I,
         ! when solving systems with salt, set up stern layer map, on the grid
         ! points

         if ( istrng /= ZERO ) call pb_ionmap( pbverbose,natom,iprob,&
              h,gox,goy,goz,xm,ym,zm,xmymzm,&
              gcrd,radi,&
              atmsas,insas,zv,iv)

         ! part II,
         ! here comes the dielectric map on the grid edges: x, y, and z

         call pb_exmol_ses( pbverbose,ifcap,ipb,natom,&
              smoothopt,dprob,epsin,epsout,&
              h,gox,goy,goz,xm,ym,zm,xmymzm,level,nfocus,&
              nwarn,nsatm,narcdot,maxarc,nbnd,nbndx,nbndy,nbndz,&
              gcrd,acrd,radi,radip3,nzratm,&
              marc,m2narc,fstarc,arcatm,dotarc,arccrd,savarc,&
              atmsas,insas,lvlset,zv,epsx,epsy,epsz,&
              iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez,savbcopt )
      end if
      call timer_stop(TIME_PBEPS)
 
      ! now call fdpb driver

      call timer_start(TIME_PBSOLV)
      call pb_fddrv( 1,natom,phi,xs(xsoffset),sbv,npbstep,npbgrid,nstlim,ionene )
      call timer_stop(TIME_PBSOLV)

      ! if requested, print a summary when the grid is set up

      if ( pbverbose .and. pbprint ) then
         call pb_print( ifcap, ipb )
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

      xsoffset = xsoffset + savxmymzm(level) + 2*savxmym(level)
!     print *, h, ionene*frcfac

   end do  !  level = 1, nfocus

   ! compute fd energy and force by the qE option
   ! note that self forces are zero
   ! delete fd coulomb energy and forces for all close pairs
   ! dbf is computed by Gilson et al

   if ( eneopt == 1 ) then

      pbfrc = ZERO ! resetting to zero for printing only
      call pb_qefrc( natom, eelrf, eelself, pbfrc, phi )

      ! when bcopt == 6, we only have reaction field energy in eelrf
      ! no need to any corrections, coulombic energy should come from
      ! pb_direct ...

      if ( savbcopt(1) /= 6 ) then 
         call pb_fdcoulomb( natom, eelcoul, pbfrc )
         eelrf = eelrf - eelself - eelcoul
      end if

!      write(6,*) 'final eelrf/ionene in kcal/mol', frcfac*eelrf, frcfac*ionene
!      write(6,*) 'final eelself in kcal/mol', frcfac*eelself
!      write(6,*) 'final eelcoul in kcal/mol', frcfac*eelcoul

      ! add ion contributions for nonlinear PBE ...

      if ( npbopt /= 0 ) then
         eelrf = eelrf + ionene
      end if

      ! return after converting to kcal/mol

      eelrf = frcfac*eelrf         

      if ( frcopt == 1 ) then
 
         ! first printing of qE forces, in electron/Angstrom^2
 
         open (unit = 102, file = 'force.dat')
         write(102,*) ' :::: Atomic qE forces ::::'
         do iatm = 1, natom
            write(102,'(3e20.6)') pbfrc(1:3,iatm)*frcfac*INV_AMBER_ELECTROSTATIC2
         end do

         pbfrc = ZERO ! resetting to zero for printing only
         call pb_dbfrc_fld(pbverbose,pbprint,natom,pbfrc,epsx,epsy,epsz,phi)

         ! second printing of DB forces, in electron/Angstrom^2

         write(102,*) ' :::: Atomic DB forces ::::'
         do iatm = 1, natom
            write(102,'(3e20.6)') pbfrc(1:3,iatm)*frcfac*INV_AMBER_ELECTROSTATIC2
         end do

      end if

   ! compute fdfrc by the charge option
   ! dbf is computed by Ye et al

   else if ( eneopt == 2 .and. epsin /= epsout ) then

      if ( frcopt == 2 ) then

         pbfrc = ZERO ! resetting to zero for printing only
         zv(1:xmymzm) = -real(insas(1:xmymzm)) ! pseudo signed distance function
         call pb_dbfrc_crg(pbverbose,pbprint,natom,eelrf,pbfrc, &
                          epsx,epsy,epsz,zv(1),phi,sbv,cphi)
      else if ( frcopt == 3 ) then
         pbfrc = ZERO ! resetting to zero for printing only
         zv(1:xmymzm) = -real(insas(1:xmymzm)) ! pseudo signed distance function
         call pb_dbfrc_fld2(pbverbose,pbprint,natom,eelrf,pbfrc, &
                           epsx,epsy,epsz,zv(1),phi,sbv,cphi)
      else
         call pb_dbene( pbverbose,pbprint,natom,eelrf,insas,phi,sbv )
      end if

   end if

 
end subroutine pb_fdfrc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Finite-difference algorithm driver
subroutine pb_fddrv( atmfirst,atmlast,phi,xs,sbv,npbstep,npbgrid,nstlim,ionene )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Solving A * x = b, where A is dielectric/salt map, b is charge map, x is phi
   ! map.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use timer_module 
   implicit none

   ! Passed variables
    
   integer atmfirst, atmlast, npbstep, npbgrid, nstlim
   _REAL_ phi(xmymzm), xs(xmymzm+2*xmym), sbv(xmymzm)
   _REAL_ ionene
    
   ! Local variables
    
   integer l
   _REAL_ factor, factor1
    
   ! initialize A and b to zero, initialize the working arrays to zero,
   ! except RD, the reciprocal of D, to ONE. save final dielectric map
   ! for boundary force calculation since the AM() array will be overwritten.
    
   bv (1:xmymzm) = ZERO
   sbv(1:xmymzm) = ZERO

   ! place charges on the grid points and save them for induced charge calculations
 
   if ( savbcopt(1) /= 6 ) then
      call pb_crggrd( atmfirst, atmlast, bv(1) )
      if ( eneopt == 2 ) then
         factor = h/frcfac*(18.2223**2)/(epsin/eps0)
         sbv(1:xmymzm) = bv(1:xmymzm)*factor
      end if
   else
      call pb_crggrd( atmfirst, atmlast, sbv(1) )
      cphi = ZERO 
      call pb_dbcgrd( bv(1), sbv(1), cphi(1), insas(1), epsx(1), epsy(1), epsz(1) )
   end if

   ! set the boundary condition, note except the first level, focusing is needed
 
   call pb_bndcnd( bv(1), sbv(1) )

   ! sbv stores the grid charge 
   if ( savbcopt(1) == 6 .and. eneopt == 2 ) then
      factor = h/frcfac*(18.2223**2)/(epsin/eps0)
      sbv(1:xmymzm) = sbv(1:xmymzm)*factor
   end if

   ! enter the core iteration routine

   call timer_start(TIME_PBITR)
   if ( npbopt == 0) then
      call solve_lpb(xm,ym,zm,xmym,xmymzm,maxitn &
                     ,fmiccg,accept,pbkappa,epsout,h,wsor &
                     ,bv(1),iv(1),xs &
                    )
   else
      call solve_npb(xm,ym,zm,xmym,xmymzm,itn,maxitn, &
                    npbstep,npbgrid, &
                    inorm,norm,wsor,lwsor, &
                    iv(1),bv(1))
   end if
   call timer_stop(TIME_PBITR)
 

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine solve_lpb(nx,ny,nz,nxny,nxnynz,p_maxitn &
                     ,p_fmiccg,p_accept,p_pbkappa,p_epsout,p_h,p_wsor &
                     ,p_bv,p_iv,p_xs &
                    )

   use pb_lsolver
   implicit none

   integer nx,ny,nz,nxny,nxnynz,p_maxitn
   _REAL_ p_fmiccg,p_accept,p_epsout,p_pbkappa,p_h,p_wsor
   _REAL_ p_bv(1:nxnynz),p_iv(1:nxnynz)
   _REAL_ p_xs(1-nxny:nxnynz+nxny)

   call init_param(nx,ny,nz,nxny,nxnynz,p_maxitn,p_fmiccg,p_accept,p_pbkappa,p_epsout,p_h,p_wsor)
   call allocate_array(solvopt)
   call init_array(solvopt,epsx,epsy,epsz,p_bv,p_iv,p_xs)

   select case ( solvopt )
   case (1)
      call pb_iccg(phi,xs)
   case (3)
      call pb_cg(phi,xs)
   case (4)
      call pb_sor(phi,xs)
   case (2)
      call pb_mg(phi,xs)
   case default
      write(6, *) 'PB bomb in pb_fddrv(): unknown solver'
      call mexit(6, 1)
   end select

   itn = l_itn
   inorm = l_inorm
   norm = l_norm

   call deallocate_array(solvopt)

end subroutine solve_lpb

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine solve_npb(nx,ny,nz,nxny,nxnynz,p_itn,p_maxitn, &
                    p_npbstep,p_npbgrid, &
                    p_inorm,p_norm,p_wsor,p_lwsor, &
                    p_iv,p_bv)

   use pb_nlsolver

   integer nx,ny,nz,nxny,nxnynz
   integer p_itn,p_maxitn,p_npbstep,p_npbgrid
   _REAL_ p_inorm,p_norm,p_wsor,p_lwsor
   _REAL_ p_iv(nxnynz),p_bv(nxnynz)

   npbstep = p_npbstep

!  if ( npbstep == 0 .or. mod(npbstep+1,npbgrid) == 0 ) then
      call init_param(nx,ny,nz,nxny,nxnynz,p_maxitn,p_npbstep,p_npbgrid, &
                      ivalence,h,pbkb,pbtemp,istrng,p_wsor,p_lwsor)

      call allocate_array(solvopt)
!  end if

   call init_array(xs,epsx,epsy,epsz,p_iv,p_bv, &
                   solvopt,npbopt,h,epsout,eps0,pbtemp,pbkappa,nbnd,iepsav )
   
   select case ( solvopt )
   case (1)
      call pb_nticcg( phi, xs, p_bv(1), accept, npbopt, nbnd, iepsav )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (2) 
      call pb_nmg( phi, xs, p_bv(1), epsout, accept, npbopt, nbnd, &
                   iepsav )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (3)
      if ( npbopt == 1 .and. ( npbstep == 0 .or. mod(npbstep+1,npbgrid) == 0 )) then
!     if ( npbopt == 1 .and. npbstep == 1 ) then
         sbv(1:xmymzm) = bv(1:xmymzm) 
         npbopt = 0
         call pb_ncg( phi, xs, laccept, npbopt )
         bv(1:xmymzm) = sbv(1:xmymzm) 
         npbopt = 1
      endif
      call pb_ncg( phi, xs, accept, npbopt )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (4)
      if ( npbopt == 1 .and. ( npbstep == 0 .or. mod(npbstep+1,npbgrid) == 0 ) ) then
!     if ( npbopt == 1 .and. npbstep == 1 ) then
         npbopt = 0
         call pb_nsor( phi, xs, lwsor, laccept, npbopt )
         npbopt = 1
      endif
      call pb_nsor( phi, xs, wsor, accept, npbopt )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (5)
      if ( npbopt == 1 .and. ( npbstep == 0 .or. mod(npbstep+1,npbgrid) == 0 ) ) then
!     if ( npbopt == 1 .and. npbstep == 1 ) then
         npbopt = 0
         call pb_asor( phi, xs, lwsor, laccept, npbopt )
         npbopt = 1
      endif
      call pb_asor( phi, xs, wsor, accept, npbopt )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (6)
      if ( npbopt == 1 .and. ( npbstep == 0 .or. mod(npbstep+1,npbgrid) == 0 ) ) then
!     if ( npbopt == 1 .and. npbstep == 1 ) then
         npbopt = 0
         call pb_dsor( phi, xs, lwsor, laccept, npbopt )
         npbopt = 1
      endif
      call pb_dsor( phi, xs, wsor, accept, npbopt )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case default
      write(6, *) 'PB bomb in pb_fddrv(): unknown solver'
      call mexit(6, 1)
   end select

   p_itn = itn
   p_inorm = inorm
   p_norm = norm

!  calculate ionic energy

   if ( npbopt /= 0 ) then
      ad (1:xmymzm) = h*ad(1:xmymzm)
      ad1(1:xmymzm) = cosh(iv(1:xmymzm)*phi(1:xmymzm)*ktinv)
      factor = h*factor/ktinv
      if ( level == nfocus ) then 
         ionene = ionene - sum(ad(1:xmymzm)*phi(1:xmymzm) + &
                              (ad1(1:xmymzm)-iv(1:xmymzm))*factor)
      else
         call pb_ionene(nfocus,level,savgox,savgoy,savgoz, &
                        savxm,savym,savzm,savh,phi,ad,ad1,iv,ionene)
      end if
   end if

!  if ( npbstep + 1 == nstlim .or. mod(npbstep+2,npbgrid) == 0 ) 
   call deallocate_array(solvopt)

end subroutine solve_npb
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
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dbcgrd( bv, sbv, cphi, insas, epsx, epsy, epsz )
   
   _REAL_ bv(xm,ym,zm)
   _REAL_ sbv(xm,ym,zm)
   _REAL_ cphi(xm,ym,zm)
   _REAL_ epsx(xm,ym,zm)
   _REAL_ epsy(xm,ym,zm)
   _REAL_ epsz(xm,ym,zm)

   integer insas(xm,ym,zm)
   integer i, j, k
   integer i0, j0, k0
   integer ip
   _REAL_ tmp, tmp0, epst
    
   tmp = ZERO; tmp0 = ZERO
   ngrdcrg = 0
   do k = 1, zm; do j = 1, ym; do i = 1, xm
      if ( sbv(i,j,k) == ZERO ) cycle
      ngrdcrg = ngrdcrg + 1

      grdcrg(1,ngrdcrg) = i
      grdcrg(2,ngrdcrg) = j
      grdcrg(3,ngrdcrg) = k
      qgrdcrg(ngrdcrg) = sbv(i,j,k)*h
   end do; end do; end do
   
   do ip = 1, nbnd
      i0 = iepsav(1,ip); j0 = iepsav(2,ip); k0 = iepsav(3,ip)
      epst = ZERO         
   
      i = i0 - 1; epst = epst + epsx(i ,j0,k0)
      if (insas(i ,j0,k0) > 0) then
         if ( cphi(i ,j0,k0) == 0.0d0 ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsx(i ,j0,k0)*cphi(i ,j0,k0)
      end if

      i = i0 + 1; epst = epst + epsx(i0,j0,k0)
      if (insas(i ,j0,k0) > 0) then
         if ( cphi(i ,j0,k0) == 0.0d0 ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsx(i0,j0,k0)*cphi(i ,j0,k0)
      end if
   
      j = j0 - 1; epst = epst + epsy(i0,j ,k0)
      if (insas(i0,j ,k0) > 0) then
         if ( cphi(i0,j ,k0) == 0.0d0 ) then
            call get_coulpot(i0,j ,k0,tmp); 
            cphi(i0,j ,k0) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsy(i0,j ,k0)*cphi(i0,j ,k0)
      end if

      j = j0 + 1; epst = epst + epsy(i0,j0,k0)
      if (insas(i0,j ,k0) > 0) then
         if ( cphi(i0,j ,k0) == 0.0d0 ) then
            call get_coulpot(i0,j ,k0,tmp)
            cphi(i0,j ,k0) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsy(i0,j0,k0)*cphi(i0,j ,k0)
      end if

      k = k0 - 1; epst = epst + epsz(i0,j0,k )
      if (insas(i0,j0,k ) > 0) then
         if ( cphi(i0,j0,k ) == 0.0d0 ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsz(i0,j0,k )*cphi(i0,j0,k )
      end if

      k = k0 + 1; epst = epst + epsz(i0,j0,k0)
      if (insas(i0,j0,k ) > 0) then
         if ( cphi(i0,j0,k ) == 0.0d0 ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsz(i0,j0,k0)*cphi(i0,j0,k )
      end if

      if (insas(i0,j0,k0) > 0) then
         if ( cphi(i0,j0,k0) == 0.0d0 ) then
            call get_coulpot(i0,j0,k0,tmp0)
            cphi(i0,j0,k0) = tmp0/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) - epst*cphi(i0,j0,k0)
      end if

      bv(i0,j0,k0) = bv(i0,j0,k0) + sbv(i0,j0,k0)
   end do
    
end subroutine pb_dbcgrd
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_coulpot(i,j,k,pot)

   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green

   integer i,j,k
   _REAL_ pot

   integer iatm
   integer itmp,jtmp,ktmp
   integer idx,idy,idz

   _REAL_ factor,qtmp,rinv

   factor = ONE/(FOURPI)/h

   pot = ZERO
   do iatm = 1, ngrdcrg
      itmp = grdcrg(1,iatm); jtmp = grdcrg(2,iatm); ktmp = grdcrg(3,iatm)
      qtmp = factor*qgrdcrg(iatm)
           
      idx = abs(i-itmp); idy = abs(j-jtmp); idz = abs(k-ktmp)
      if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
         rinv = green(idx,idy,idz)
         pot = pot + qtmp*rinv
      else
         rinv = ONE/sqrt(REAL(idx**2 + idy**2 + idz**2))
         pot = pot + qtmp*rinv
      end if
   end do  !  iatm = 1, ngrdcrg

end subroutine get_coulpot
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign Debye-Huckel potential for the boundary charge grid
subroutine pb_bndcnd( bv, sbv )
   
   ! Common variables
    
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
     
   ! Passed variables
    
   _REAL_ bv(xm,ym,zm), sbv(xm,ym,zm)
    
   ! Local variables
    
   integer i, j, k, iatm
   integer xmtmp, ymtmp, zmtmp, ix, iy, iz, itmp, jtmp, ktmp, idx, idy, idz
   _REAL_ htmp, goxtmp, goytmp, goztmp
   _REAL_ qtmp
   _REAL_ x, y, z
   _REAL_ xi, yi, zi, aa, bb, cc, aa1, bb1, cc1
   _REAL_ r, rinv
   _REAL_ factor
    
   factor = ONE/(FOURPI)
!  ix = 0; iy = 0; iz = 0
    
   ! bcopt = 1
   ! use zero potential for boundary grid, the boundary will be all solvent
    
   if ( bcopt == 1 ) then
   !  write(6, *) "PB bomb in pb_bndcnd(): zero BC not supported"
   !  call mexit(6, 1)
    
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

   ! bcopt = 6
   ! sum of grid charge reaction field contribution. the boundary will be all solvent.
    
   else if ( bcopt == 6 ) then
       
      ngrdcrg = 0
      do k = 1, zm; do j = 1, ym; do i = 1, xm
         if ( sbv(i,j,k) == ZERO ) cycle
         ngrdcrg = ngrdcrg + 1
         grdcrg(1,ngrdcrg) = i
         grdcrg(2,ngrdcrg) = j
         grdcrg(3,ngrdcrg) = k
         qgrdcrg(ngrdcrg) = sbv(i,j,k)*h
      end do; end do; end do
       
      factor = factor/h !* epsout * (ONE/epsout - ONE/epsin) 

      do iatm = 1, ngrdcrg
         itmp = grdcrg(1,iatm); jtmp = grdcrg(2,iatm); ktmp = grdcrg(3,iatm)
         qtmp = factor*qgrdcrg(iatm)
           
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
 
         !if (i==1 .and. j==1) write(6,*) idx,idy,idz,pbkappa,qtmp,rinv
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


end subroutine pb_fddrv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ total finite difference es energy and forces
subroutine pb_qefrc( natom, grdnrg, grdself, pbfrc, phi )
   
   implicit none

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
   
   implicit none

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
              gci1*( gcj1*l_green(dijx ,dijy ,dijz ) + gcj2*l_green(dijx1,dijy ,dijz ) + &
                     gcj3*l_green(dijx ,dijy1,dijz ) + gcj4*l_green(dijx1,dijy1,dijz ) + &
                     gcj5*l_green(dijx ,dijy ,dijz1) + gcj6*l_green(dijx1,dijy ,dijz1) + &
                     gcj7*l_green(dijx ,dijy1,dijz1) + gcj8*l_green(dijx1,dijy1,dijz1) ) + &
              gci2*( gcj1*l_green(dijx2,dijy ,dijz ) + gcj2*l_green(dijx ,dijy ,dijz ) + &
                     gcj3*l_green(dijx2,dijy1,dijz ) + gcj4*l_green(dijx ,dijy1,dijz ) + &
                     gcj5*l_green(dijx2,dijy ,dijz1) + gcj6*l_green(dijx ,dijy ,dijz1) + &
                     gcj7*l_green(dijx2,dijy1,dijz1) + gcj8*l_green(dijx ,dijy1,dijz1) ) + &
              gci3*( gcj1*l_green(dijx ,dijy2,dijz ) + gcj2*l_green(dijx1,dijy2,dijz ) + &
                     gcj3*l_green(dijx ,dijy ,dijz ) + gcj4*l_green(dijx1,dijy ,dijz ) + &
                     gcj5*l_green(dijx ,dijy2,dijz1) + gcj6*l_green(dijx1,dijy2,dijz1) + &
                     gcj7*l_green(dijx ,dijy ,dijz1) + gcj8*l_green(dijx1,dijy ,dijz1) ) + &
              gci4*( gcj1*l_green(dijx2,dijy2,dijz ) + gcj2*l_green(dijx ,dijy2,dijz ) + &
                     gcj3*l_green(dijx2,dijy ,dijz ) + gcj4*l_green(dijx ,dijy ,dijz ) + &
                     gcj5*l_green(dijx2,dijy2,dijz1) + gcj6*l_green(dijx ,dijy2,dijz1) + &
                     gcj7*l_green(dijx2,dijy ,dijz1) + gcj8*l_green(dijx ,dijy ,dijz1) ) + &
              gci5*( gcj1*l_green(dijx ,dijy ,dijz2) + gcj2*l_green(dijx1,dijy ,dijz2) + &
                     gcj3*l_green(dijx ,dijy1,dijz2) + gcj4*l_green(dijx1,dijy1,dijz2) + &
                     gcj5*l_green(dijx ,dijy ,dijz ) + gcj6*l_green(dijx1,dijy ,dijz ) + &
                     gcj7*l_green(dijx ,dijy1,dijz ) + gcj8*l_green(dijx1,dijy1,dijz ) ) + &
              gci6*( gcj1*l_green(dijx2,dijy ,dijz2) + gcj2*l_green(dijx ,dijy ,dijz2) + &
                     gcj3*l_green(dijx2,dijy1,dijz2) + gcj4*l_green(dijx ,dijy1,dijz2) + &
                     gcj5*l_green(dijx2,dijy ,dijz ) + gcj6*l_green(dijx ,dijy ,dijz ) + &
                     gcj7*l_green(dijx2,dijy1,dijz ) + gcj8*l_green(dijx ,dijy1,dijz ) ) + &
              gci7*( gcj1*l_green(dijx ,dijy2,dijz2) + gcj2*l_green(dijx1,dijy2,dijz2) + &
                     gcj3*l_green(dijx ,dijy ,dijz2) + gcj4*l_green(dijx1,dijy ,dijz2) + &
                     gcj5*l_green(dijx ,dijy2,dijz ) + gcj6*l_green(dijx1,dijy2,dijz ) + &
                     gcj7*l_green(dijx ,dijy ,dijz ) + gcj8*l_green(dijx1,dijy ,dijz ) ) + &
              gci8*( gcj1*l_green(dijx2,dijy2,dijz2) + gcj2*l_green(dijx ,dijy2,dijz2) + &
                     gcj3*l_green(dijx2,dijy ,dijz2) + gcj4*l_green(dijx ,dijy ,dijz2) + &
                     gcj5*l_green(dijx2,dijy2,dijz ) + gcj6*l_green(dijx ,dijy2,dijz ) + &
                     gcj7*l_green(dijx2,dijy ,dijz ) + gcj8*l_green(dijx ,dijy ,dijz ) ) )
        
         ffx = dble( gcij( 1)*(l_green(dijx1,dijy ,dijz ) - l_green(dijx2,dijy ,dijz )) + &
                     gcij( 2)*(l_green(dijx0,dijy ,dijz ) - l_green(dijx ,dijy ,dijz )) + &
                     gcij( 3)*(l_green(dijx1,dijy1,dijz ) - l_green(dijx2,dijy1,dijz )) + &
                     gcij( 4)*(l_green(dijx0,dijy1,dijz ) - l_green(dijx ,dijy1,dijz )) + &
                     gcij( 5)*(l_green(dijx1,dijy ,dijz1) - l_green(dijx2,dijy ,dijz1)) + &
                     gcij( 6)*(l_green(dijx0,dijy ,dijz1) - l_green(dijx ,dijy ,dijz1)) + &
                     gcij( 7)*(l_green(dijx1,dijy1,dijz1) - l_green(dijx2,dijy1,dijz1)) + &
                     gcij( 8)*(l_green(dijx0,dijy1,dijz1) - l_green(dijx ,dijy1,dijz1)) + &
                     gcij( 9)*(l_green(dijx ,dijy ,dijz ) - l_green(dijx3,dijy ,dijz )) + &
                     gcij(10)*(l_green(dijx ,dijy1,dijz ) - l_green(dijx3,dijy1,dijz )) + &
                     gcij(11)*(l_green(dijx ,dijy ,dijz1) - l_green(dijx3,dijy ,dijz1)) + &
                     gcij(12)*(l_green(dijx ,dijy1,dijz1) - l_green(dijx3,dijy1,dijz1)) + &
                     gcij(13)*(l_green(dijx1,dijy2,dijz ) - l_green(dijx2,dijy2,dijz )) + &
                     gcij(14)*(l_green(dijx0,dijy2,dijz ) - l_green(dijx ,dijy2,dijz )) + &
                     gcij(15)*(l_green(dijx1,dijy2,dijz1) - l_green(dijx2,dijy2,dijz1)) + &
                     gcij(16)*(l_green(dijx0,dijy2,dijz1) - l_green(dijx ,dijy2,dijz1)) + &
                     gcij(17)*(l_green(dijx ,dijy2,dijz ) - l_green(dijx3,dijy2,dijz )) + &
                     gcij(18)*(l_green(dijx ,dijy2,dijz1) - l_green(dijx3,dijy2,dijz1)) + &
                     gcij(19)*(l_green(dijx1,dijy ,dijz2) - l_green(dijx2,dijy ,dijz2)) + &
                     gcij(20)*(l_green(dijx0,dijy ,dijz2) - l_green(dijx ,dijy ,dijz2)) + &
                     gcij(21)*(l_green(dijx1,dijy1,dijz2) - l_green(dijx2,dijy1,dijz2)) + &
                     gcij(22)*(l_green(dijx0,dijy1,dijz2) - l_green(dijx ,dijy1,dijz2)) + &
                     gcij(23)*(l_green(dijx ,dijy ,dijz2) - l_green(dijx3,dijy ,dijz2)) + &
                     gcij(24)*(l_green(dijx ,dijy1,dijz2) - l_green(dijx3,dijy1,dijz2)) + &
                     gcij(25)*(l_green(dijx1,dijy2,dijz2) - l_green(dijx2,dijy2,dijz2)) + &
                     gcij(26)*(l_green(dijx0,dijy2,dijz2) - l_green(dijx ,dijy2,dijz2)) + &
                     gcij(27)*(l_green(dijx ,dijy2,dijz2) - l_green(dijx3,dijy2,dijz2)) )
       
         ffy = dble( gcij( 1)*(l_green(dijx ,dijy1,dijz ) - l_green(dijx ,dijy2,dijz )) + &
                     gcij( 2)*(l_green(dijx1,dijy1,dijz ) - l_green(dijx1,dijy2,dijz )) + &
                     gcij( 3)*(l_green(dijx ,dijy0,dijz ) - l_green(dijx ,dijy ,dijz )) + &
                     gcij( 4)*(l_green(dijx1,dijy0,dijz ) - l_green(dijx1,dijy ,dijz )) + &
                     gcij( 5)*(l_green(dijx ,dijy1,dijz1) - l_green(dijx ,dijy2,dijz1)) + &
                     gcij( 6)*(l_green(dijx1,dijy1,dijz1) - l_green(dijx1,dijy2,dijz1)) + &
                     gcij( 7)*(l_green(dijx ,dijy0,dijz1) - l_green(dijx ,dijy ,dijz1)) + &
                     gcij( 8)*(l_green(dijx1,dijy0,dijz1) - l_green(dijx1,dijy ,dijz1)) + &
                     gcij( 9)*(l_green(dijx2,dijy1,dijz ) - l_green(dijx2,dijy2,dijz )) + &
                     gcij(10)*(l_green(dijx2,dijy0,dijz ) - l_green(dijx2,dijy ,dijz )) + &
                     gcij(11)*(l_green(dijx2,dijy1,dijz1) - l_green(dijx2,dijy2,dijz1)) + &
                     gcij(12)*(l_green(dijx2,dijy0,dijz1) - l_green(dijx2,dijy ,dijz1)) + &
                     gcij(13)*(l_green(dijx ,dijy ,dijz ) - l_green(dijx ,dijy3,dijz )) + &
                     gcij(14)*(l_green(dijx1,dijy ,dijz ) - l_green(dijx1,dijy3,dijz )) + &
                     gcij(15)*(l_green(dijx ,dijy ,dijz1) - l_green(dijx ,dijy3,dijz1)) + &
                     gcij(16)*(l_green(dijx1,dijy ,dijz1) - l_green(dijx1,dijy3,dijz1)) + &
                     gcij(17)*(l_green(dijx2,dijy ,dijz ) - l_green(dijx2,dijy3,dijz )) + &
                     gcij(18)*(l_green(dijx2,dijy ,dijz1) - l_green(dijx2,dijy3,dijz1)) + &
                     gcij(19)*(l_green(dijx ,dijy1,dijz2) - l_green(dijx ,dijy2,dijz2)) + &
                     gcij(20)*(l_green(dijx1,dijy1,dijz2) - l_green(dijx1,dijy2,dijz2)) + &
                     gcij(21)*(l_green(dijx ,dijy0,dijz2) - l_green(dijx ,dijy ,dijz2)) + &
                     gcij(22)*(l_green(dijx1,dijy0,dijz2) - l_green(dijx1,dijy ,dijz2)) + &
                     gcij(23)*(l_green(dijx2,dijy1,dijz2) - l_green(dijx2,dijy2,dijz2)) + &
                     gcij(24)*(l_green(dijx2,dijy0,dijz2) - l_green(dijx2,dijy ,dijz2)) + &
                     gcij(25)*(l_green(dijx ,dijy ,dijz2) - l_green(dijx ,dijy3,dijz2)) + &
                     gcij(26)*(l_green(dijx1,dijy ,dijz2) - l_green(dijx1,dijy3,dijz2)) + &
                     gcij(27)*(l_green(dijx2,dijy ,dijz2) - l_green(dijx2,dijy3,dijz2)) )
      
         ffz = dble( gcij( 1)*(l_green(dijx ,dijy ,dijz1) - l_green(dijx ,dijy ,dijz2)) + &
                     gcij( 2)*(l_green(dijx1,dijy ,dijz1) - l_green(dijx1,dijy ,dijz2)) + &
                     gcij( 3)*(l_green(dijx ,dijy1,dijz1) - l_green(dijx ,dijy1,dijz2)) + &
                     gcij( 4)*(l_green(dijx1,dijy1,dijz1) - l_green(dijx1,dijy1,dijz2)) + &
                     gcij( 5)*(l_green(dijx ,dijy ,dijz0) - l_green(dijx ,dijy ,dijz )) + &
                     gcij( 6)*(l_green(dijx1,dijy ,dijz0) - l_green(dijx1,dijy ,dijz )) + &
                     gcij( 7)*(l_green(dijx ,dijy1,dijz0) - l_green(dijx ,dijy1,dijz )) + &
                     gcij( 8)*(l_green(dijx1,dijy1,dijz0) - l_green(dijx1,dijy1,dijz )) + &
                     gcij( 9)*(l_green(dijx2,dijy ,dijz1) - l_green(dijx2,dijy ,dijz2)) + &
                     gcij(10)*(l_green(dijx2,dijy1,dijz1) - l_green(dijx2,dijy1,dijz2)) + &
                     gcij(11)*(l_green(dijx2,dijy ,dijz0) - l_green(dijx2,dijy ,dijz )) + &
                     gcij(12)*(l_green(dijx2,dijy1,dijz0) - l_green(dijx2,dijy1,dijz )) + &
                     gcij(13)*(l_green(dijx ,dijy2,dijz1) - l_green(dijx ,dijy2,dijz2)) + &
                     gcij(14)*(l_green(dijx1,dijy2,dijz1) - l_green(dijx1,dijy2,dijz2)) + &
                     gcij(15)*(l_green(dijx ,dijy2,dijz0) - l_green(dijx ,dijy2,dijz )) + &
                     gcij(16)*(l_green(dijx1,dijy2,dijz0) - l_green(dijx1,dijy2,dijz )) + &
                     gcij(17)*(l_green(dijx2,dijy2,dijz1) - l_green(dijx2,dijy2,dijz2)) + &
                     gcij(18)*(l_green(dijx2,dijy2,dijz0) - l_green(dijx2,dijy2,dijz )) + &
                     gcij(19)*(l_green(dijx ,dijy ,dijz ) - l_green(dijx ,dijy ,dijz3)) + &
                     gcij(20)*(l_green(dijx1,dijy ,dijz ) - l_green(dijx1,dijy ,dijz3)) + &
                     gcij(21)*(l_green(dijx ,dijy1,dijz ) - l_green(dijx ,dijy1,dijz3)) + &
                     gcij(22)*(l_green(dijx1,dijy1,dijz ) - l_green(dijx1,dijy1,dijz3)) + &
                     gcij(23)*(l_green(dijx2,dijy ,dijz ) - l_green(dijx2,dijy ,dijz3)) + &
                     gcij(24)*(l_green(dijx2,dijy1,dijz ) - l_green(dijx2,dijy1,dijz3)) + &
                     gcij(25)*(l_green(dijx ,dijy2,dijz ) - l_green(dijx ,dijy2,dijz3)) + &
                     gcij(26)*(l_green(dijx1,dijy2,dijz ) - l_green(dijx1,dijy2,dijz3)) + &
                     gcij(27)*(l_green(dijx2,dijy2,dijz ) - l_green(dijx2,dijy2,dijz3)) )
      
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

contains

function l_green (i,j,k)

   implicit none
   integer i,j,k
   _REAL_ l_green

   if ( i <= 20  .and. j <= 20 .and. k <= 20 ) then
      l_green = green(i,j,k)
   else
      l_green = ONE/sqrt(dble(i*i+j*j+k*k))
   end if

end function l_green
   
end subroutine pb_fdcoulomb 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary energy and forces
subroutine pb_dbene( verbose,eneout,natom,eel,insas,phi,sbv )

   use solvent_accessibility, only : dprob, radi, arccrd, dotarc, arcatm
   implicit none
 
   ! Passed variables
 
   logical verbose, eneout
   integer natom
   integer insas(xm,ym,zm)
   _REAL_ phi(xm,ym,zm), sbv(xm,ym,zm)
   _REAL_ eel
 
   ! Local variables
 
   integer  ip, i, j, k, iatm, jatm
   _REAL_ srfcrg, factor, scalfact
   _REAL_ g(3), x(3), dx(3), crd(3), dx2, dist, rdist, acg
   _REAL_ dinv, de, eelrf, pcolomb

   _REAL_, parameter :: smallcrg = 0.5d0
 
   pcolomb = ZERO

   ! calculate the total potential in solute if bcopt=6
   if ( savbcopt(1) == 6 ) then
      do k = 1, zm
         do j = 1, ym
            do i = 1, xm
               if ( insas(i,j,k) > 0 ) then
                  call get_coulpot(i,j,k,pcolomb)
                  phi(i,j,k) = phi(i,j,k) + pcolomb / eps0
               end if
            end do
         end do
      end do
   end if
   
   ! initialization
 
   factor = THREE*h/(TWOPI)
   scalfact = ONE
 
   srfcrg = ZERO; eel = ZERO

   ! for InsightII display
   !open (unit=55, file='ms.dot')
   !write (55, '("DOTS")')
 
   do ip = 1, nbnd

      ! collecting boundary grid info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k

      ! project the boundary grid point on to the molecular surface

      if      ( abs(insas(i,j,k)) == 2 .and. iatm > 0 ) then
         ! the contact boundary grid points are projected to the atom spheres
         x(1:3) = acrd(1:3,iatm)
         dist = radi(iatm)
      else if ( abs(insas(i,j,k)) == 1 .and. iatm > 0) then
         ! the reentry boundary grid points are projected to the solvent spheres
         x(1:3) = arccrd(1:3,iatm)
         dist = dprob
      else
         ! otherwise do not project. The setting should make crd = g
         x(1:3) = g(1:3)
         dist = ONE
      end if

      dx = g - x; dx2 = dx(1)**2 + dx(2)**2 + dx(3)**2
      if ( dx2 == ZERO ) then
         rdist = ONE
      else
         rdist = dist*ONE/sqrt(dx2)
      end if
      crd = x + dx*rdist

      ! for InsightII display
      !write (55,'(4(f8.3,2x))') crd(1:3), 300.
 
      ! compute induced charge on the molecular surface

      acg = factor*(phi(i,j,k)-&
      SIXTH*( phi(i-1,j,k)+phi(i+1,j,k)+phi(i,j-1,k)+phi(i,j+1,k)+phi(i,j,k-1)+phi(i,j,k+1) ))
      acg = acg - sbv(i,j,k)
      srfcrg = srfcrg + acg*frcfac/(AMBER_ELECTROSTATIC2)
 
      ! compute reaction field energy due to this boundary grid point

      eelrf = ZERO
      do jatm = 1, natom
         dx(1) = crd(1) - acrd(1,jatm)
         dx(2) = crd(2) - acrd(2,jatm)
         dx(3) = crd(3) - acrd(3,jatm)
         dinv = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)

         de = acg*acrg(jatm)*dinv
         eelrf = eelrf + de
      end do
      eel = eel + eelrf
 
   end do ! end of ip = 1, nbnd

   ! for InsightII display
   !close(55)
   !stop

   if ( scalerf .and. abs(totcrg) > smallcrg ) then
      scalfact = abs( totcrg/srfcrg*(ONE/epsin - ONE/epsout)*eps0 )
      srfcrg = scalfact*srfcrg
      eel = scalfact*eel
   end if
    
   eel = HALF*frcfac*eel
    
   !if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eel
   !end if

contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_coulpot(i,j,k,pot)

   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green

   integer i,j,k
   _REAL_ pot

   integer iatm
   integer itmp,jtmp,ktmp
   integer idx,idy,idz

   _REAL_ factor,qtmp,rinv

   factor = ONE/(FOURPI)/h

   pot = ZERO
   do iatm = 1, ngrdcrg
      itmp = grdcrg(1,iatm); jtmp = grdcrg(2,iatm); ktmp = grdcrg(3,iatm)
      qtmp = factor*qgrdcrg(iatm)
           
      idx = abs(i-itmp); idy = abs(j-jtmp); idz = abs(k-ktmp)
      if (idx <= 20 .and. idy <= 20 .and. idz <= 20) then
         rinv = green(idx,idy,idz)
         pot = pot + qtmp*rinv
      else
         rinv = ONE/sqrt(REAL(idx**2 + idy**2 + idz**2))
         pot = pot + qtmp*rinv
      end if
   end do  !  iatm = 1, ngrdcrg

end subroutine get_coulpot

end subroutine pb_dbene
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dbfrc_fld(verbose,eneout,natom,f,epsx,epsy,epsz,phi)
  
   use solvent_accessibility, only : radi, arccrd, dotarc, arcatm
   implicit none

#  include "constants.h"

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ f(3,natom)
   _REAL_ epsx(xm,ym,zm), epsy(xm,ym,zm), epsz(xm,ym,zm), phi(xm,ym,zm)

   ! local variables

   integer i, j, k, iatm, matm, natm, iarc
   integer xp, yp, zp
   _REAL_ df, factor, factor1, sfactor
   _REAL_ lEx, lEy, lEz
   _REAL_ dx, dy, dz 
   _REAL_ dfx, dfy, dfz 
   _REAL_ adx, ady, adz
   _REAL_ x(3), crd(3)
   _REAL_ d1, d2 
    
   ! initialization

   sfactor = -HALF*(epsout-epsin)/(epsin*epsout)

   ! contributions from epsx boundary edges

   do xp = 1, nbndx
      i = iepsavx(1,xp); j = iepsavx(2,xp); k = iepsavx(3,xp); iatm = iepsavx(4,xp)
      lEx = phi(i+1,j,k) - phi(i,j,k)
      crd(1) = gox + h*(i+fedgex(xp)); crd(2) = goy + h*j; crd(3) = goz + h*k
      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
      else
         x(1:3) = arccrd(1:3,-iatm)
      end if
      dx = crd(1) - x(1)
      dy = crd(2) - x(2)
      dz = crd(3) - x(3)
      adx = abs(dx)

      df = sfactor*lEx*lEx*epsx(i,j,k)**2
      factor1 = df/adx

#     include "pb_dbfrc_fld.h"

   enddo

   ! contributions from epsy boundary grid edges
             
   do yp = 1, nbndy
      i = iepsavy(1,yp); j = iepsavy(2,yp); k = iepsavy(3,yp); iatm = iepsavy(4,yp)
      lEy = phi(i,j+1,k) - phi(i,j,k)
      crd(1) = gox + h*i; crd(2) = goy + h*(j+fedgey(yp)); crd(3) = goz + h*k

      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
      else
         x(1:3) = arccrd(1:3,-iatm)
      end if

      dx = crd(1) - x(1)
      dy = crd(2) - x(2)
      dz = crd(3) - x(3)
      ady = abs(dy)

      df = sfactor*lEy*lEy*epsy(i,j,k)**2
      factor1= df/ady

#     include "pb_dbfrc_fld.h"

   enddo
           
   ! contributions from epsz boundary grid edges
           
   do zp = 1, nbndz
      i = iepsavz(1,zp); j = iepsavz(2,zp); k = iepsavz(3,zp); iatm = iepsavz(4,zp)
      lEz = phi(i,j,k+1) - phi(i,j,k)
      crd(1) = gox + h*i; crd(2) = goy + h*j; crd(3) = goz + h*(k+fedgez(zp))
      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
      else
         x(1:3) = arccrd(1:3,-iatm)
      end if

      dx = crd(1) - x(1)
      dy = crd(2) - x(2)
      dz = crd(3) - x(3)
      adz = abs(dz)

      df = sfactor*lEz*lEz*epsz(i,j,k)**2
      factor1= df/adz

#     include "pb_dbfrc_fld.h"

   end do

end subroutine pb_dbfrc_fld
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dbfrc_crg( verbose,eneout,natom,eelrf,f,epsx,epsy,epsz,insas,phi,sbv,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd, dotarc, arcatm
   implicit none

#  include "constants.h"

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ eelrf
   _REAL_ f(3,natom)
   _REAL_ epsx(xm,ym,zm), epsy(xm,ym,zm), epsz(xm,ym,zm)
   _REAL_ insas(xm,ym,zm), phi(xm,ym,zm), sbv(xm,ym,zm), cphi(xm,ym,zm)

   ! local variables

   integer, parameter :: n_point = 7
   integer i, j, k, iatm, jatm, matm, natm, iarc, ip
   _REAL_ srfcrg
   _REAL_ g(3), x(3), dx(3), crd(3), dist, rdist, crg
   _REAL_ dum(3), coulomb(3)
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ fx0, fy0, fz0
   _REAL_ dudxi0, dudyi0, dudzi0
   _REAL_ d1, d2
   _REAL_, allocatable :: pol_charge(:)

   allocate (pol_charge(1:nbnd))

   ! initialization

   !mjhsieh added for silencing gfortran
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO; x = ZERO; dist = ZERO
   coulomb = ZERO; srfcrg = ZERO

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg

   call get_charge_pol(nbnd,phi,cphi,sbv,pol_charge,srfcrg)

   ! main double loops over polarization charges and atom charges

   do ip = 1, nbnd

      ! collecting boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
      crg = pol_charge(ip)

      ! project the surface grid point on to the molecular surface, crd() is the
      ! new coord, and x() is the atom/probe coord, fx/y/z0 is the grid version
      ! of crd()

      if ( iatm == 0 ) then
         write(6,*) 'PB Bomb in pb_dbfrc_fld(): can not find projection atom/probe'
         call mexit(6, 1)
      end if

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

      ! grid unit of the projected boundary point

      fx0 = (crd(1)-gox)/h
      fy0 = (crd(2)-goy)/h
      fz0 = (crd(3)-goz)/h

      ! inner loop over atoms...
      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and natom atomic charges

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      ! now it is time to compute dbf
      ! compute E on the inner side of surface position crd(1:3)

      call gradu(xm,ym,zm,-ONE,n_point,4,fx0,fy0,fz0,dudxi0,dudyi0,dudzi0,phi,cphi,insas)

      dudxi0 = -dudxi0*FOURPI*eps0
      dudyi0 = -dudyi0*FOURPI*eps0
      dudzi0 = -dudzi0*FOURPI*eps0
    
      ! add the coulomb field to get the total E of inner side  

      dudxi0 = dudxi0 + coulomb(1)
      dudyi0 = dudyi0 + coulomb(2)
      dudzi0 = dudzi0 + coulomb(3)

      ! compute the DBF as HALF*Q*D

      dum(1) = HALF*crg*dudxi0
      dum(2) = HALF*crg*dudyi0
      dum(3) = HALF*crg*dudzi0
      
      ! collecting contact forces

      if ( abs(insas(i,j,k)) == 2 .and. iatm > 0 ) then
         f(1,iatm) = f(1,iatm) - dum(1)
         f(2,iatm) = f(2,iatm) - dum(2)
         f(3,iatm) = f(3,iatm) - dum(3)

      ! collecting reentry forces

      else if ( abs(insas(i,j,k)) == 1 .and. iatm > 0 ) then
         iarc = dotarc(iatm)
         matm = arcatm(1,iarc); natm = arcatm(2,iarc)
         d1=abs(sqrt((crd(1)-acrd(1,matm))**2+(crd(2)-acrd(2,matm))**2+(crd(3)-acrd(3,matm))**2)-radi(matm))
         d2=abs(sqrt((crd(1)-acrd(1,natm))**2+(crd(2)-acrd(2,natm))**2+(crd(3)-acrd(3,natm))**2)-radi(natm))
         f(1:3,matm) = f(1:3,matm) - d2/(d1+d2)*dum(1:3) 
         f(1:3,natm) = f(1:3,natm) - d1/(d1+d2)*dum(1:3) 
      end if

   end do

   deallocate (pol_charge)

   open (unit = 103, file = 'force.dat')

   write(103,*) ' :::: Atomic qE forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') qex(iatm),qey(iatm),qez(iatm)
   end do

   write(103,*) ' :::: Atomic DB forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') f(1:3,iatm)
   end do

   eelrf = eelrf * AMBER_ELECTROSTATIC2 / 2.d0
   !if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf!*AMBER_ELECTROSTATIC2
   !end if

end subroutine pb_dbfrc_crg
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary energy and forces
subroutine pb_dbfrc_fld2( verbose,eneout,natom,eel,f,epsx,epsy,epsz,insas,phi,sbv,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd, dotarc, arcatm

#  include "constants.h"   

   ! Passed variables

   logical verbose, eneout
   integer natom
   _REAL_ insas(xm,ym,zm), phi(xm,ym,zm), sbv(xm,ym,zm), cphi(xm,ym,zm)
   _REAL_ epsx(xm,ym,zm), epsy(xm,ym,zm), epsz(xm,ym,zm)
   _REAL_ eel, f(3,natom)

   ! Local variables

   integer  i, j, k, iatm, matm, natm, iarc, ip
   _REAL_ x(3), crd(3), dum(3)
   _REAL_ crg, eelrf
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ fx(natom), fy(natom), fz(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ fx0, fy0, fz0
   _REAL_ dudxi0, dudyi0, dudzi0
   _REAL_ dudxo0, dudyo0, dudzo0
   _REAL_ rn(1:3), dn
   _REAL_ dudni, dudno, coulomb(3) 
   _REAL_ dr, ds, total_s, f_bnd
   _REAL_ EN, rsphere, d1, d2

   ! initialization for DBF

   ax = acrd(1,1:natom)
   ay = acrd(2,1:natom)
   az = acrd(3,1:natom)

   fx = ZERO; fy = ZERO; fz = ZERO

   do ip = 1, nbndx
      i = iepsavx(1,ip); j = iepsavx(2,ip); k = iepsavx(3,ip); iatm = iepsavx(4,ip)
      fx0 = i + fedgex(ip); fy0 = j; fz0 = k
      crd(1) = gox + h*i + fedgex(ip)*h; crd(2) = goy + h*j; crd(3) = goz + h*k

      if (iatm == 0) then
         write(6,*) 'PBMD FATAL ERROR: can not find projection atom/probe' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
      else 
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
      end if

      if ( iatm > 0) then
         rn(1) = crd(1) - x(1)
         rn(2) = crd(2) - x(2)
         rn(3) = crd(3) - x(3)
      else
         rn(1) =  x(1) - crd(1)
         rn(2) =  x(2) - crd(2)
         rn(3) =  x(3) - crd(3)
      end if
      dr = abs(rn(1))

#     include "pb_dbfrc_fld2.h"

   end do !nbndx

   do ip = 1, nbndy
      i = iepsavy(1,ip); j = iepsavy(2,ip); k = iepsavy(3,ip); iatm = iepsavy(4,ip)
      fx0 = i; fy0 = j + fedgey(ip); fz0 = k
      crd(1) = gox + h*i; crd(2) = goy + h*j + fedgey(ip)*h; crd(3) = goz + h*k

      if (iatm == 0) then
         write(6,*) 'PBMD FATAL ERROR: can not find projection atom/probe' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
      else 
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
      end if

      if ( iatm > 0) then
         rn(1) = crd(1) - x(1)
         rn(2) = crd(2) - x(2)
         rn(3) = crd(3) - x(3)
      else
         rn(1) =  x(1) - crd(1)
         rn(2) =  x(2) - crd(2)
         rn(3) =  x(3) - crd(3)
      end if
      dr = abs(rn(2))

#     include "pb_dbfrc_fld2.h"

   end do !nbndy

   do ip = 1, nbndz
      i = iepsavz(1,ip); j = iepsavz(2,ip); k = iepsavz(3,ip); iatm = iepsavz(4,ip)
      fx0 = i; fy0 = j; fz0 = k + fedgez(ip)
      crd(1) = gox + h*i ; crd(2) = goy + h*j; crd(3) = goz + h*k+fedgez(ip)*h

      if (iatm == 0) then
         write(6,*) 'PBMD FATAL ERROR: can not find projection atom/probe' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
      else 
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
      end if

      if ( iatm > 0 ) then
         rn(1) = crd(1) - x(1)
         rn(2) = crd(2) - x(2)
         rn(3) = crd(3) - x(3)
      else
         rn(1) =  x(1) - crd(1)
         rn(2) =  x(2) - crd(2)
         rn(3) =  x(3) - crd(3)
      end if
      dr = abs(rn(3))

#     include "pb_dbfrc_fld2.h"

   end do !nbndz

   open (unit = 103, file = 'force.dat')

   write(103,*) ' :::: Atomic qE forces ::::'
   qex = ZERO; qey = ZERO; qez = ZERO
   do iatm = 1, natom
      write(103,'(3e20.6)') qex(iatm),qey(iatm),qez(iatm)
   end do

   write(103,*) ' :::: Atomic DB forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') f(1:3,iatm)
   end do

   !if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total molecular surface', total_s
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf!*AMBER_ELECTROSTATIC2
   !end if

end subroutine pb_dbfrc_fld2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine gradu(l,m,n,phip,n2,n3,xp,yp,zp,dudx,dudy,dudz,u,cu,phi)

   ! passed variables

   integer l,m,n
   integer n2,n3
   _REAL_  xp,yp,zp,phip
   _REAL_  u(1:l,1:m,1:n), cu(1:l,1:m,1:n), phi(1:l,1:m,1:n)
   _REAL_  dudx,dudy,dudz

   ! local variable

   integer job,info
   integer ix0,iy0,iz0,ix1,iy1,iz1,ix,iy,iz
   integer nsub
   integer i1, i2, j2, k2, i, j, k
   integer, parameter :: nn = 100
   _REAL_  w1(1:n3,1:nn)
   !_REAL_  w2(1:n3,1:nn)
   _REAL_  b(1:nn)
   _REAL_  sd(1:n3+1)
   _REAL_  ew(1:nn),w4(1:nn)
   _REAL_  w3ux(1:nn),w3uy(1:nn),w3uz(1:nn)
   _REAL_  ewux(1:n3),ewuy(1:n3),ewuz(1:n3)
   _REAL_  uw(1:n3,1:n3)
   _REAL_  vl(1:nn,1:nn)
   !_REAL_  uxw(1:n3,1:n2)
   !_REAL_  uxwxv(1:n2,1:n2)
   _REAL_  dist, dist0
   _REAL_  dx,dy,dz

   !mjhsieh added for silencing gfortran
   ix = 0; iy = 0; iz = 0
   !w1(1:n3,1:nn) = 0.0d0; w2(1:n3,1:nn) = 0.0d0

   ! select the grid (ix,iy,iz) closest to the surface charge (xp,yp,zp)
   
   ix0 = nint(xp)
   iy0 = nint(yp)
   iz0 = nint(zp)
   dist0 = 9999.0d0
   do i = -1 ,1
      do j = -1 ,1
         do k = -1 ,1
            ix1 = ix0 + i
            iy1 = iy0 + j
            iz1 = iz0 + k
            if ( phi(ix1,iy1,iz1) > ZERO ) cycle
!           dist = (fx0-ix1)**2+(fy0-iy1)**2+(fz0-iz1)**2
            dist = (xp-ix1)**2+(yp-iy1)**2+(zp-iz1)**2
            if ( dist < dist0 ) then
               ix = ix1
               iy = iy1
               iz = iz1
               dist0 = dist
            end if
         end do
      end do
   end do

   ! select the closest inside grid points to interplate E_in  
   ! and construcut the matrix for SVD

   nsub = 0
   i1 = 0
   do while ( nsub < n2 )
      do i2 = -i1, i1
         do j2 = -i1, i1
            do k2 = -i1, i1
               i = i2 + ix
               j = j2 + iy
               k = k2 + iz
               if ( i > l .or. i < 1 ) cycle
               if ( j > m .or. j < 1 ) cycle
               if ( k > n .or. k < 1 ) cycle
               if ( i2*i2 + j2*j2 + k2*k2 == i1 ) then
                  if ( phi(i,j,k)*phip >= ZERO ) then
                     nsub = nsub + 1
                     dx = (i - xp)*h
                     dy = (j - yp)*h
                     dz = (k - zp)*h
                     if ( nsub <= n2 ) then
                        w1(1,nsub) = ONE
                        w1(2,nsub) = dx
                        w1(3,nsub) = dy
                        w1(4,nsub) = dz
                        b(nsub) = u(i,j,k)
                     end if       
                 end if
              end if
            end do
         end do
      end do
      i1 = i1 + 1
   end do
   if ( nsub > n2 ) nsub = n2

!  do i = 1, nsub
!     write(120,'(5x,5e20.6)') w1(1:4, i), b(i)
!  end do
!  do i =1,n2
!     w3ux(i) = 0.0
!     w3uy(i) = 0.0
!     w3uz(i) = 0.0
!  end do
!  do i =1,n3
!     ewux(i) = 0.0
!     ewuy(i) = 0.0
!     ewuz(i) = 0.0
!  end do

   ! call dsvdc for the singular value decomposition

!  print *,'n2=',n2
!  print *,'n3=',n3
!  print *,'b'
!  print *,v
!  w2 = w1
   job = 11

   call dsvdc(w1,n3,n3,nsub,sd,ew,uw,n3,vl,nn,w4,job,info)

!  print *,v
!  print *,'SVD is completed'

!  print *,'~~~~~~~~~~~~ u ~~~~~~~~~~~'
!  do i = 1, n3
!     do j = 1, n3
!        print *,sum(uw(1:n3,i)*uw(1:n3,j))
!     end do
!  end do
!  print *,'~~~~~~~~~~~~ u ~~~~~~~~~~~'

!  print *,'~~~~~~~~~~~~ v ~~~~~~~~~~~'
!  do i = 1, n2
!     do j = 1, n2
!        print *,sum(vl(1:n2,i)*vl(1:n2,j))
!     end do
!  end do
!  print *,'~~~~~~~~~~~~ v ~~~~~~~~~~~'

!  do i = 1, n3
!     do j = 1, n2
!        uxw(i,j)=sum(uw(1:n3,i)*w2(1:n3,j))
!     end do
!  end do
!  print *,'~~~~~~~~~~~~~~~~~~~~~~~~~'
!  do i = 1, n3
!     do j = 1, n2
!        uxwxv(i,j) = sum(uxw(i,1:n2)*vl(1:n2,j))
!     end do
!     print "(100f10.3)",uxwxv(i,1:n2)
!  end do
!  print *,'~~~~~~~~~~~~~~~~~~~~~~~~~'

   ! calculate E_in using the returned least-squared coefficients 

   do i = 1, n3-1
      ewux(i)= uw(2,i)/sd(i)
      ewuy(i)= uw(3,i)/sd(i)
      ewuz(i)= uw(4,i)/sd(i)
   enddo

   if ( sd(n3) > 1.0d-12 ) then
      ewux(n3) = uw(2,n3)/sd(n3)
      ewuy(n3) = uw(3,n3)/sd(n3)
      ewuz(n3) = uw(4,n3)/sd(n3)
   else
      ewux(n3) = ZERO
      ewuy(n3) = ZERO
      ewuz(n3) = ZERO
   endif
!  print *,v

   do i = 1, nsub
      w3ux(i) = ZERO
      w3uy(i) = ZERO
      w3uz(i) = ZERO
      do j = 1, n3
         w3ux(i) = w3ux(i) + vl(i,j)*ewux(j)
         w3uy(i) = w3uy(i) + vl(i,j)*ewuy(j)
         w3uz(i) = w3uz(i) + vl(i,j)*ewuz(j)
      enddo
!     print *,i,v
   enddo
!  print *,v

   dudx = ZERO
   dudy = ZERO
   dudz = ZERO
   do i = 1, nsub
      dudx = dudx +  w3ux(i)*b(i)
      dudy = dudy +  w3uy(i)*b(i)
      dudz = dudz +  w3uz(i)*b(i)
   end do
!  print *,v


end subroutine gradu
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_charge_pol(nbnd,phi,cphi,sbv,pol_charge,srfcrg)

#  include "constants.h"

   ! passed variables
   
   integer nbnd  
   _REAL_ phi(xm,ym,zm), cphi(xm,ym,zm), sbv(xm,ym,zm)
   _REAL_ pol_charge(nbnd), srfcrg

   ! local variables

   integer i, j, k, ip
   _REAL_ total_phi(7), factor, total

   factor = THREE*h/(TWOPI)

   total = ZERO
   do ip = 1, nbnd

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip)

      total_phi(1) = phi(i,  j,k) + cphi(i,  j,k)  
      total_phi(2) = phi(i-1,j,k) + cphi(i-1,j,k)  
      total_phi(3) = phi(i+1,j,k) + cphi(i+1,j,k)  
      total_phi(4) = phi(i,j-1,k) + cphi(i,j-1,k)  
      total_phi(5) = phi(i,j+1,k) + cphi(i,j+1,k)  
      total_phi(6) = phi(i,j,k-1) + cphi(i,j,k-1)  
      total_phi(7) = phi(i,j,k+1) + cphi(i,j,k+1)  

      pol_charge(ip) = factor*( total_phi(1)-&
           SIXTH*( total_phi(2)+total_phi(3)+&
                   total_phi(4)+total_phi(5)+&
                   total_phi(6)+total_phi(7) ) )
      pol_charge(ip) = ( pol_charge(ip) - sbv(i,j,k))*frcfac*INV_AMBER_ELECTROSTATIC2
      total = total + pol_charge(ip) 
   end do
   srfcrg = total


end subroutine get_charge_pol

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_coulomb(natom,crg,crd,ac,ax,ay,az,fx,fy,fz,eelrf,coulomb)

   ! passed variables

   integer natom
   _REAL_ crg, crd(1:3), coulomb(1:3), eelrf
   _REAL_ ac(natom), ax(natom), ay(natom), az(natom)
   _REAL_ fx(natom), fy(natom), fz(natom)

   ! local variables

   integer jatm
   _REAL_ dinv, d2inv, de, dx(1:3), dff, dcc

   coulomb(1:3) = ZERO
   do jatm = 1, natom
      dx(1) = crd(1) - ax(jatm)
      dx(2) = crd(2) - ay(jatm)
      dx(3) = crd(3) - az(jatm)
      dinv = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2); d2inv = dinv**2

      de = ac(jatm)*dinv
      dcc = de*d2inv
      dff = crg*dcc

      ! calculate reaction field energy

      eelrf = eelrf + de*crg

      ! calculate the QE force on the atom

      fx(jatm) = fx(jatm) + dx(1)*dff
      fy(jatm) = fy(jatm) + dx(2)*dff
      fz(jatm) = fz(jatm) + dx(3)*dff

      ! calculate the coulomb field on the surface charge

      coulomb(1) = coulomb(1) + dx(1)*dcc
      coulomb(2) = coulomb(2) + dx(2)*dcc
      coulomb(3) = coulomb(3) + dx(3)*dcc
   end do


end subroutine get_coulomb


end module poisson_boltzmann
