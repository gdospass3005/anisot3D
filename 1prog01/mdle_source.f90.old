MODULE mdle_source

  IMPLICIT NONE

CONTAINS

  SUBROUTINE source_dervgauss_init(trs,dt,ns,fpeak,tdelay)
    IMPLICIT NONE
    REAL, INTENT(INOUT), DIMENSION(ns) :: trs
    INTEGER, INTENT(IN) :: ns
    REAL, INTENT(IN) :: dt, fpeak
    REAL, INTENT(OUT) :: tdelay
    INTEGER :: i
    REAL :: bw,t,pi=3.1415926536

    tdelay = 4.0 /(3.*fpeak*SQRT(2.*pi))

    bw=1.
    do i=1,ns
       t=(i-1)*dt - tdelay
       trs(i) = -dervgauss(t,fpeak,bw)
    end do

    OPEN(81,FILE='sourcememory.ad',STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)
    WRITE(81, REC=1) trs(:)
    CLOSE(81)

  END SUBROUTINE source_dervgauss_init


  REAL FUNCTION dervgauss(t,fpeak,bombweight)
    REAL, INTENT(IN) :: t,fpeak
    REAL :: x,xx,pi=3.1415926536
    REAL, INTENT(IN) :: bombweight

    x=3*fpeak*t
    xx=x*x

    dervgauss= bombweight*(-2.*pi*t)*EXP(-pi*xx)

  END FUNCTION dervgauss


END MODULE mdle_source

