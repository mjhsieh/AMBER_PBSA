#include "copyright.h"
#  define _REAL_ double precision
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ The PBSA main program
program pbsamain
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  Meng-Juei Hsieh
   !  The Luo Research Group
   !  University of California, Irvine
   !
   !  Setup MPI and file handling. 
   !  Call pbsa to perform calculations.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

#  include "files.h"

   numgroup = 1
   call pbsafile()
   
   call pbsa()

   call mexit(6,0)

end program pbsamain

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ The PBSA driver
subroutine pbsa()
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  Meng-Juei Hsieh
   !  The Luo Research Group
   !  University of California, Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use timer_module
   use decomp, only : allocate_int_decomp, allocate_real_decomp, &
                      deallocate_int_decomp, deallocate_real_decomp

   implicit none

   logical belly, erstop
   integer ier,ncalls

#  include "files.h"
#  include "memory.h"
#  include "box.h"
#  include "md.h"
#  include "parms.h"
#  include "extra.h"
#  include "timer.h"

   _REAL_ ene(51)
   integer nr
   _REAL_ carrms
   character(len=8) initial_date, setup_end_date, final_date
   character(len=10) initial_time, setup_end_time, final_time
   _REAL_,  allocatable :: x(:)
   integer, allocatable :: ix(:), ipairs(:)
   character(len=4), allocatable :: ih(:)
   _REAL_  r_stack(1)
   integer i_stack(1)

   ier = 0
   
   ! Initialize the cpu timer. Needed for machines where returned cpu times
   ! are relative.

   call date_and_time( initial_date, initial_time )
   call timer_init()
    
   ! in the single-threaded version, the one process is master

   master = .true.

   erstop = .false.

   ! Only the master node (only node when single-process)
   ! performs the initial setup and reading/writing

   call timer_start(TIME_TOTAL)
   if ( master ) then

      !        --- first, initial reads to determine memry sizes:

      call timer_start(TIME_READ)
      call mdread1()
!     call openparm(parm//CHAR(0))
!     call rdparm1()
      call myopen(8,parm,'O','F','R')
      call rdparm1(8)

      !        --- now, we can allocate memory:

      call locmem

      !        --- dynamic memory allocation:

      REQUIRE(lastr>0.and.lasti>0.and.lastpr>0.and.lasth>0)
      allocate( x      (lastr), stat = ier ); REQUIRE(ier==0)
      allocate( ix     (lasti), stat = ier ); REQUIRE(ier==0)
      allocate( ipairs(lastpr), stat = ier ); REQUIRE(ier==0)
      allocate( ih     (lasth), stat = ier ); REQUIRE(ier==0)

      if( idecomp > 0 ) then
         call allocate_int_decomp(natom, nres)
      else
         call allocate_int_decomp(1, 1)
      endif

      lastrst = 1
      lastist = 1
      r_stack(1) = 0.0d0
      i_stack(1) = 0

      write(6,'(/,a,5x,a)') '|','Memory Use     Allocated'
      write(6,'(a,5x,a,i14)') '|', 'Real      ', lastr
      write(6,'(a,5x,a,i14)') '|', 'Hollerith ', lasth
      write(6,'(a,5x,a,i14)') '|', 'Integer   ', lasti
      write(6,'(a,5x,a,i14)') '|', 'Max Pairs ', lastpr
      write(6,'(a,5x,a,i14)') '|', 'Max Rstack', lastrst
      write(6,'(a,5x,a,i14)') '|', 'Max Istack', lastist
      write(6,'(a,5x,a,i14,a)') '|', '  Total   ', &
           (8*(lastr+lastrst) + 4*(lasth+lasti+lastpr+lastist))/1024, ' kbytes'

      !        --- second reads to finalize

!     call rdparm2(x,ix,ih)
      call rdparm2(x,ix,ih,ipairs,8,i_stack)
      call mdread2(x,ix,ih)

!        --- alloc memory for decomp module that needs info from mdre
!        ad2
      if( idecomp == 1 .or. idecomp == 2 ) then
         call allocate_real_decomp(nres)
      else if( idecomp == 3 .or. idecomp == 4 ) then
         call allocate_real_decomp(npdec*npdec)
      end if

      ! EVALUATE SOME CONSTANTS FROM MDREAD SETTINGS

      nr = nrp
      belly = ibelly > 0

      ! READ COORDINATES AND VELOCITIES

      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,t)
      if ( belly ) call bellyf()
      call timer_stop(TIME_READ)
      
      ! Compatibility print-out please ignore

      write(6,'(" Number of triangulated 3-point waters found: ",i8)') 0

      ! OPEN THE DATA DUMPING FILES AND POSITION IT DEPENDING ON THE TYPE OF RUN

      ! call open_dump_files

      if ( master ) call amflsh(6)

      ! end of master process setup
   end if  ! (master)
   
   call date_and_time( setup_end_date, setup_end_time )

   ! Now do the dynamics or minimization.

   if ( master ) write(6,'(/80(1H-)/''   4.  RESULTS'',/80(1H-)/)')

   ! Input flag imin determines the type of calculation: MD, minimization, ...

   select case ( imin )
   case ( 0 )
      ! Dynamics:

      call timer_start(TIME_RUNMD)
      call runmd(x,ix,ih,ipairs, &
              x(lcrd),x(lwinv),x(lmass),x(lforce), &
              x(lvel),x(lvel2),x(l45),x(lcrdr), &
              x(l50),x(l95),ix(i70),x(l75),erstop,r_stack,i_stack)
      call timer_stop(TIME_RUNMD)

      if (master) call amflsh(6)

   case ( 1 )
      ! Minimization:

      ! Input flag ntmin determines the method of minimization
      select case ( ntmin )
      case ( 0, 1, 2 )
         call runmin(x,ix,ih,ipairs,x(lcrd),x(lforce),x(lvel), &
                 ix(iibh),ix(ijbh),x(l50),x(lwinv),ix(lasti), &
                 x(l95),ene,r_stack,i_stack, carrms)
      case default
         ! invalid ntmin
         ! ntmin input validation occurs in mdread.f
         write(6,'(/2x,a,i3,a)') 'Error: Invalid NTMIN (',ntmin,').'
         call mexit(6,0)
      end select

      if (master) call amflsh(6)

   case ( 5 )
      ! trajene option not supported
      ! imin input validation is in mdread.f
      write (6,*) 'Error: Post-processing of trajectory not supported'
      call mexit(6,0)
   case default
      ! invalid imin
      ! imin input validation is in mdread.f
      write(6,'(/2x,a,i3,a)') 'Error: Invalid IMIN (',imin,').'
      call mexit(6,0)
   end select

   !     -- calc time spent running vs setup

   call timer_stop(TIME_TOTAL)
   call date_and_time( final_date, final_time )
   call timer_summary()

   if ( master ) then

      !     --- write out final times

      write(6,'(12(a))') '|           Job began  at ', initial_time(1:2), &
           ':', initial_time(3:4), ':', initial_time(5:10), '  on ',&
           initial_date(5:6), '/', initial_date(7:8), '/', initial_date(1:4)
      write(6,'(12(a))') '|           Setup done at ', setup_end_time(1:2),  &
           ':', setup_end_time(3:4), ':', setup_end_time(5:10), '  on ', &
           setup_end_date(5:6), '/',setup_end_date(7:8),'/',setup_end_date(1:4)
      write(6,'(12(a))') '|           Run   done at ', final_time(1:2),  &
           ':', final_time(3:4), ':', final_time(5:10), '  on ', &
           final_date(5:6), '/', final_date(7:8), '/', final_date(1:4)
      call nwallclock( ncalls )
      write(6, '(''|'',5x,''wallclock() was called'',I8,'' times'')') ncalls
   end if

   !     --- dynamic memory deallocation:

   call pb_free()
   deallocate( ih, stat = ier ); REQUIRE(ier==0)
   deallocate( ipairs, stat = ier ); REQUIRE(ier==0)
   deallocate( ix, stat = ier ); REQUIRE(ier==0)
   deallocate( x, stat = ier ); REQUIRE(ier==0)


end subroutine pbsa

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ The PBSA file handler
subroutine pbsafile
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  Meng-Juei Hsieh
   !  The Luo Research Group
   !  University of California, Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

#  include "files.h"

   character(len=80) argbuf
   integer iarg ! index of the current argument
   integer pb_iargc ! wrapper to intrinsic that returns the index of the last argument
                    ! from either the command line or a string
   integer last_arg_index ! index of the last argument

   !     --- default file names ---

   mdin   = 'mdin'
   mdout  = 'mdout'
   inpcrd = 'inpcrd'
   parm   = 'prmtop'
   restrt = 'restrt'
   refc   = 'refc'
   mdvel  = 'mdvel'
   mden   = 'mden'
   mdcrd  = 'mdcrd'
   mdinfo = 'mdinfo'
   vecs   = 'vecs'
   freqe  = 'dummy'
   rstdip = 'rstdip'
   inpdip = 'inpdip'
   mddip  = 'mddip'
   radii  = 'radii'
   cpin   = 'cpin'
   cpout  = 'cpout'
   cprestrt = 'cprestrt'
   if (numgroup == 1) groups = ' '

   !     --- default status of output: New

   owrite = 'N'

   !     --- get command line arguments ---

   iarg = 0
   last_arg_index = pb_iargc()
   do while (iarg < last_arg_index)
      iarg = iarg + 1

      call getpb_arg(iarg,argbuf)
      select case (argbuf)
      case (' ')
         continue
      case ('-O')
#ifdef ABSOFT_WINDOWS
         owrite = 'U' !      status of output: unknown
#else
         owrite = 'R' !      status of output: Replace
#endif
      case ('-i')
         iarg = iarg + 1
         call getpb_arg(iarg,mdin)
      case ('-o')
         iarg = iarg + 1
         call getpb_arg(iarg,mdout)
      case ('-p')
         iarg = iarg + 1
         call getpb_arg(iarg,parm)
      case ('-c')
         iarg = iarg + 1
         call getpb_arg(iarg,inpcrd)
      case ('-vec')
         iarg = iarg + 1
         call getpb_arg(iarg,vecs)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-radii')
         iarg = iarg + 1
         call getpb_arg(iarg,radii)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-f')
         iarg = iarg + 1
         call getpb_arg(iarg,freqe)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-r')
         iarg = iarg + 1
         call getpb_arg(iarg,restrt)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-ref', '-z')
         iarg = iarg + 1
         call getpb_arg(iarg,refc)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-e')
         iarg = iarg + 1
         call getpb_arg(iarg,mden)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-v')
         iarg = iarg + 1
         call getpb_arg(iarg,mdvel)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-x', '-t')
         iarg = iarg + 1
         call getpb_arg(iarg,mdcrd)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-inf')
         iarg = iarg + 1
         call getpb_arg(iarg,mdinfo)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-idip')
         iarg = iarg + 1
         call getpb_arg(iarg,inpdip)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-rdip')
         iarg = iarg + 1
         call getpb_arg(iarg,rstdip)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-mdip')
         iarg = iarg + 1
         call getpb_arg(iarg,mddip)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-cpin')
         iarg = iarg + 1
         call getpb_arg(iarg,cpin)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-cpout')
         iarg = iarg + 1
         call getpb_arg(iarg,cpout)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-cprestrt')
         iarg = iarg + 1
         call getpb_arg(iarg,cprestrt)
         write(6,*) "Unsupported command arguments, exit.";call mexit(6,0)
      case default
         write(6,'(/,5x,a,a)') 'ERROR: Unknown file: ',argbuf
         write(6,9000)
         call mexit(6, 1)
      end select 
   end do  !  while (iarg < last_arg_index)
 

   9000 format(/,5x, &
         'usage: pbsa  [-O] -i mdin -o mdout -p prmtop -c inpcrd ', &
         '-r restrt',/19x,'[-ref refc -x mdcrd -v mdvel -e mden ', &
         '-idip inpdip -rdip rstdip -mdip mddip ', &
         '-inf mdinfo -radii radii]' &
         , /, 'Consult the manual for additional options.')

end subroutine pbsafile

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ RESERVED FOR FUTURE MPI EXTENSION
integer function pb_iargc()
   implicit none
   integer iargc
   pb_iargc = iargc()

end function pb_iargc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ RESERVED FOR FUTURE MPI EXTENSION
subroutine getpb_arg(iarg, arg)
   implicit none
   integer iarg
   character(len=*) arg
   ! Intrinsic getarg requires a 4 byte integer argument; 
   ! this guards the argument for builds with default 8 byte integers.
   call getarg(int(iarg,4), arg)
 
end subroutine getpb_arg

!  Rewritten by: Meng-Juei Hsieh
subroutine mexit(filenum, exitstatus)
   implicit none
   integer filenum
   integer exitstatus

   if (filenum > 0) then  ! close this unit if greater than zero
      close(unit=filenum)
   endif
   ! exit status; error if non-zero
#if XLF90 || IBM3090 || F2C
   if (exitstatus.ne.0) stop 1; stop 0
#else
   call exit(exitstatus)
#endif
end subroutine mexit
