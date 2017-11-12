MODULE mdle_source

IMPLICIT NONE

CONTAINS
SUBROUTINE source_init(trs,dt,ns,fpeak,tdelay)
  IMPLICIT NONE
  REAL, INTENT(INOUT), DIMENSION(ns) :: trs
  INTEGER, INTENT(IN) :: ns
  REAL, INTENT(IN) :: dt, fpeak
  REAL, INTENT(OUT) :: tdelay
  INTEGER :: i
  REAL :: t,pi=3.141592653589793238462643383279502884197

  REAL :: wpeak,waux,tt

  wpeak = 2.*pi*fpeak
  waux  = 0.5*wpeak

   tdelay = 6./(5.*fpeak)


  do i=1,ns
     t=(i-1)*dt
     tt = t - tdelay

     trs(i) = exp(-waux*waux*tt*tt/4.)*cos(wpeak*tt)

  end do

END SUBROUTINE source_init

END MODULE
