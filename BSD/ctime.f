#define _REAL_ double precision
!#######################################################################

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ctime here]
subroutine ctime(time)
   !#######################################################################
   !     RETURNS CURRENT CPU-TIME IN SECONDS (REAL*4!)
   
   !      REAL ACCUM
   !      INTEGER          IERROR, ITICKS, MCLOCK
   !      REAL             CPUTIM
   
   !MPQ  CALL TIMER(ITICKS)
   !MPQ  TIME=0.01*REAL(ITICKS)
   
   !===> IBM/VS-Fortran
   
   !      CALL CPUTIME (ACCUM, IERROR)
   !      TIME = ACCUM * 1.0D-6
   
   !===> IBM-AIX
   
   !      TIME = timer() * 0.01
   
   !===> Ardent/Titan
   
   !ARD  TIME = CPUTIM (0.0)
   
   !===> Targon / iPSC/2
   
   !TARG TIME = MCLOCK () * 1E-3
   
   return
end subroutine ctime 
