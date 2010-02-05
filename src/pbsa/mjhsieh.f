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

subroutine prepb_init(mycn1,mycn2,mynttyp)
   implicit none
#  include "parms.h"
   _REAL_ mycn1(*),mycn2(*)
   integer mynttyp
   cn1(1:mynttyp)=mycn1(1:mynttyp)
   cn2(1:mynttyp)=mycn2(1:mynttyp)
   return
end subroutine prepb_init

subroutine mypb_force(natom,nres,ntypes,npdec,ipres,iac,ico,exclat,&
                   cn1,cn2,cg,xx,f,epol,evdw,eelt,esurf,edisp)
   use poisson_boltzmann, only : pb_force
   use dispersion_cavity, only : np_force
   use timer_module
   implicit none
#  include "md.h"
#  include "timer.h"

   integer natom,   nres,  ntypes,npdec
   integer ipres(*),iac(*),ico(*),exclat(*)
   _REAL_  cn1(*),  cn2(*),cg(*), xx(*),    f(*)
   _REAL_  evdw,eelt,epol,esurf,edisp
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
   if ( inp /= 0 ) then
      esurf = 0.0d0; edisp = 0.0d0
      call np_force(natom,nres,ntypes,ipres,iac,ico, &
                    cn1,cn2,xx,f,esurf,edisp)
   end if
   call timer_stop(TIME_TOTAL)
   call date_and_time( final_date, final_time )
   return
end subroutine mypb_force
!ene(2) = evdw
!ene(3) = eelt
!ene(4) = epol
