! <compile=optimized>
#include "copyright.h"
#  define _REAL_ double precision

subroutine mexit(filenum, exitstatus)
   implicit none
   integer filenum
   integer exitstatus

   if (filenum > 0) then  ! close this unit if greater than zero
      close(unit=filenum)
   endif
   ! exit status; error if non-zero
#  if XLF90 || IBM3090 || F2C
   if (exitstatus.ne.0) stop 1; stop 0
#  else
   call exit(exitstatus)
#  endif
end subroutine mexit

subroutine private_getx(natom,filenum,x)
   implicit none
   integer natom, filenum
   _REAL_ x(*)

   integer nr3, ier, i
   character(len=80) title1
   character(len=256) testbuffer
   nr3=3*natom
   call myopen(filenum,'unittest.inpcrd','O','F','R')
   read(filenum,'(a80)') title1
   read(filenum,'(a80)') testbuffer
   if( testbuffer(6:6) == ' ' ) then ! this is an old file
      read(testbuffer,'(i5,e15.7)') natom
   else                              ! assume a new file
      read(testbuffer,'(i6,e15.7)') natom
   end if
   read(filenum,'(6f12.7)') (x(i),i=1,nr3)
   close(filenum, iostat=ier)
   return
end subroutine private_getx

subroutine mypb_force(natom,nres,ntypes,ipres,iac,ico,exclat,&
                   cn1,cn2,cg,xx,f,epol)
   use poisson_boltzmann, only : pb_force
   use timer_module
   implicit none
#  include "md.h"
#  include "timer.h"

   integer natom,   nres,  ntypes
   integer ipres(*),iac(*),ico(*),exclat(*)
   _REAL_  cn1(*),  cn2(*),cg(*), xx(*),    f(*)
   _REAL_  evdw,eelt,epol
   character(len=8) initial_date, setup_end_date, final_date
   character(len=10) initial_time, setup_end_time, final_time
   ! Initialize the cpu timer. Needed for machines where returned cpu times
   call date_and_time( initial_date, initial_time )
   call timer_init()
   call timer_start(TIME_TOTAL)
   if( ipb /= 0 ) then
      call pb_force(natom,nres,ntypes,ipres,iac,ico,exclat, &
                    cn1,cn2,cg,xx,f,evdw,eelt,epol)
   end if
   call timer_stop(TIME_TOTAL)
   call date_and_time( final_date, final_time )
   return
end subroutine mypb_force
!ene(2) = evdw
!ene(3) = eelt
!ene(4) = epol

#ifdef NIL
subroutine mjhsieh(xx,ix,x,f,ener,vir)


   integer ix(*)

#  include "pb_constants.h"
#  include "md.h"
#  include "memory.h"
#  include "parms.h"
#  include "pb_md.h"


   _REAL_  enmr(3)

   _REAL_  x(*),f(*),ene(30),vir(*)
   _REAL_  ener(*)

   integer i,m
   _REAL_  evdw,eelt,e3bod,epol,esurf
   _REAL_  epolar,aveper,aveind,avetot


   ! ZERO OUT THE ENERGIES AND FORCES

   aveper=0.d0
   aveind=0.d0
   avetot=0.d0
   !dipiter=0.d0
   !dvdl=0.d0
   !dipole_temp=0.d0
   do i=1,3
      enmr(i) = 0.d0
   end do
   do i=1,4
      vir(i) = 0.d0
   end do
   !virvsene = 0.d0
   do i=1,3*natom
      f(i) = 0.d0
   end do

   epolar = 0.d0
   e3bod = 0.d0

   ! part I: bonded terms
   ! part II: restraining terms
   ! part III: implicit solvent nonbonded treatments

   esurf = 0.d0

   ! pb options


   ! part IV: summary of energy components for printing
   !
   !    ene(1):    total energy
   !    ene(2):    van der Waals
   !    ene(3):    electrostatic energy
   !    ene(4):    10-12 (hb) energy, or GB/PB energy when igb.gt.0
   !    ene(5):    bond energy
   !    ene(6):    angle energy
   !    ene(7):    torsion angle energy
   !    ene(8):    1-4 nonbonds
   !    ene(9):    1-4 electrostatics
   !    ene(10):   constraint energy
   !    ene(11-19):  used a scratch, but not needed further below
   !    ene(20):   position constraint energy
   !    ene(21):   charging free energy result
   !    ene(22):   noe volume penalty
   !    ene(23):   surface-area dependent solvation energy or cavity energy
   !    ene(24):   surface-area dependent dispersion energy

   do m = 2,15
      ene(1) = ene(1) + ene(m)
   end do
   ene(1) = ene(1) + epolar + e3bod + ene(23) + ene(24)

   ene(5) = ene(6)+ene(7)
   ene(6) = ene(8)+ene(9)
   ene(7) = ene(10)+ene(13)
   ene(8) = ene(11)+ene(14)
   ene(9) = ene(12)+ene(15)
   ene(10) = ene(17)+ene(20)+enmr(1)+enmr(2)+enmr(3)
   ene(1) = ene(1)+ene(10)

   ! transfer the energy array to external usage for printing

   ener(1:10) = ene(1:10)
   ener(11) = epolar
   ener(12) = aveper
   ener(13) = aveind
   ener(14) = avetot
   ener(15) = ene(23)
   ener(16) = e3bod
   ener(17) = ene(21)
   ener(18) = ene(24)

end subroutine mjhsieh
#endif
