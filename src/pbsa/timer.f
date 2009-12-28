#include "copyright.h"
#include "../include/dprec.fh"
#include "timer.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ SYSTEM CLOCK
subroutine wallclock( wallc )
   implicit none
   _REAL_ wallc
   integer ncalls,n
   integer mycount, myrate

   data ncalls /0/

   call system_clock( COUNT=mycount, COUNT_RATE=myrate)
   wallc = dble(mycount)/dble(myrate)
   ncalls = ncalls + 1
   return
entry nwallclock ( n )
   n = ncalls
   return
end subroutine wallclock

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBSA TIMER
module timer_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Meng-Juei Hsieh
!  The Luo Research Group
!  University of California, Irvine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none
   logical timerswitch(maxtime)
   _REAL_ time_bulletin (maxtime)
   character(len=25) time_desc (maxtime)

contains

subroutine timer_init
   implicit none
   timerswitch=.false.
   time_bulletin=0d0

   time_desc(TIME_TOTAL)   ='Total time'
   time_desc(TIME_READ)    ='Read prm/crd time'
   time_desc(TIME_RUNMD)   ='Runmd Time'
   time_desc(TIME_FORCE)   ='Force time'
   time_desc(TIME_BOND)    ='Bond/Angle/Dihedral'
   time_desc(TIME_NONBON)  ='Nonbond force'
   !-- pb
   time_desc(TIME_PBFORCE) ='PB Nonbond'
   time_desc(TIME_PBLIST)  ='PB NB list'
   time_desc(TIME_PBSETUP) ='PB FD grid'
   time_desc(TIME_PBSAS)   ='PB Sasa'
   time_desc(TIME_PBFDFRC) ='PB FD Force'
   time_desc(TIME_PBEPS)   ='PB Epsmap'
   time_desc(TIME_PBSOLV)  ='PB Solver'
   time_desc(TIME_PBITR)   ='PB ITER'
   time_desc(TIME_PBDBE)   ='PB DB Force'
   time_desc(TIME_PBDIRECT)='PB Direct'
   time_desc(TIME_PBMP)    ='PB Multiple'
   time_desc(TIME_SINH)    ='PB SINH'
   !-- np
   time_desc(TIME_NPFORCE) ='NP Nonbond'
   time_desc(TIME_NPSAS)   ='NP Sasa'
   time_desc(TIME_NPCAV)   ='NP Cavity'
   time_desc(TIME_NPDIS)   ='NP Dispersion'
   return
end subroutine timer_init

subroutine timer_start( num_event )
   implicit none

   integer num_event
   _REAL_ mytime

   call wallclock(mytime)
   if ( num_event > maxtime ) then
      write(6,*)'index ',num_event,' bigger than MAXTIME ',maxtime
      write(6,*)'attempt to add timer '
      call mexit(6,1)
   else if ( timerswitch(num_event) ) then
      write(6,*)'timer: need to reset the time before start'
      call mexit(6,1)
   else if (time_bulletin(num_event) <= 0) then
      time_bulletin(num_event) = mytime
   else
      time_bulletin(num_event) = mytime-time_bulletin(num_event)
   endif
   timerswitch(num_event) = .true.
   return
end subroutine timer_start

subroutine timer_stop( num_event )
   implicit none
   integer num_event
   _REAL_ mytime

   call wallclock(mytime)
   if ( num_event > maxtime )then
      write(6,*)'index ',num_event,' bigger than MAXTIME ',maxtime
      write(6,*)'attempt to close timer '
      call mexit(6,1)
   else if ( .not. timerswitch(num_event) ) then
      write(6,*)'timer: need to start the time before stop'
      call mexit(6,1)
   end if
   time_bulletin(num_event) = mytime-time_bulletin(num_event)
   timerswitch(num_event) = .false.
   return
end subroutine timer_stop

subroutine timer_summary
   implicit none
   integer i

   write(6,'(/80(1H-)/,''   5.  TIMINGS'',/80(1H-)/)')
   do i = 1, maxtime
      if ( .not. timerswitch(i) ) then
         if (time_bulletin(i) > 0) &
            write(6,"('|',1x,a,f10.2)") time_desc(i),time_bulletin(i)
      else
         write(6,*) "Warning, this timer is not stopped: ",time_desc(i)
      endif
   enddo
   return
end subroutine timer_summary


end module timer_module
