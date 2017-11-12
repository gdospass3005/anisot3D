PROGRAM model3D

  USE mdle_model3D
  IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! VARIABLES USED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER :: i,j,Nx,Ny,Nz,Io
  REAL :: d1,d2
  REAL :: da1,da2
  REAL :: db1,db2
  REAL :: dc1,dc2
  REAL :: dd1,dd2
  REAL :: de1,de2

  REAL, DIMENSION(:,:,:), ALLOCATABLE :: panel

  CHARACTER(LEN=20) :: infile='data.dat'
  CHARACTER(LEN=20) :: outfile1='C11.ad'
  CHARACTER(LEN=20) :: outfile2='C33.ad'
  CHARACTER(LEN=20) :: outfile3='C44.ad'
  CHARACTER(LEN=20) :: outfile4='C13.ad'
  CHARACTER(LEN=20) :: outfile5='rhox.ad'
  CHARACTER(LEN=20) :: outfile6='C66.ad'

  LOGICAL :: verbose = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  OPEN(25,FILE=infile,STATUS='UNKNOWN')

  READ(25,*) Nx, Ny, Nz
  READ(25,*) Io
  READ(25,*) d1,d2
  READ(25,*) da1,da2
  READ(25,*) db1,db2
  READ(25,*) dc1,dc2
  READ(25,*) dd1,dd2
  READ(25,*) de1,de2

  PRINT*,'PROGRAMA GERADOR DE PAINEIS PARA A MODELAGEM 3-D EM MEIOS TI'

  IF (verbose .EQV. .TRUE.) THEN
     PRINT*, Nx,Ny,Nz,' = grid size in x, y and z'
     PRINT*, Io,' = depth of the plane interface in grid units'
     PRINT*, d1,d2,' = C11'
     PRINT*, da1,da2,' = C33'
     PRINT*, db1,db2,' = C44'
     PRINT*, dc1,dc2,' = C13'
     PRINT*, dd1,dd2,' = rhox'
     PRINT*, de1,de2,' = C66'
  END IF

  ALLOCATE(panel(Nx,Ny,Nz))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL model3Dsub(Io,panel,d1,d2)
  OPEN(35,FILE=outfile1,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO J=1,Ny  
  DO I=1,Nx  
     WRITE(35, REC=i+(j-1)*Nx) panel(i,j,:)
  END DO
  END DO
  CLOSE(35)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C33
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL model3Dsub(Io,panel,da1,da2)
  OPEN(35,FILE=outfile2,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO J=1,Ny  
  DO I=1,Nx  
     WRITE(35, REC=i+(j-1)*Nx) panel(i,j,:)
  END DO
  END DO
  CLOSE(35)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C44
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL model3Dsub(Io,panel,db1,db2)
  OPEN(35,FILE=outfile3,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO J=1,Ny  
  DO I=1,Nx  
     WRITE(35, REC=i+(j-1)*Nx) panel(i,j,:)
  END DO
  END DO
  CLOSE(35)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C13
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL model3Dsub(Io,panel,dc1,dc2)
  OPEN(35,FILE=outfile4,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO J=1,Ny  
  DO I=1,Nx  
     WRITE(35, REC=i+(j-1)*Nx) panel(i,j,:)
  END DO
  END DO
  CLOSE(35)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! rhox
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL model3Dsub(Io,panel,dd1,dd2)
  OPEN(35,FILE=outfile5,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO J=1,Ny 
  DO I=1,Nx  
     WRITE(35, REC=i+(j-1)*Nx) panel(i,j,:)
  END DO
  END DO
  CLOSE(35)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! C66
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL model3Dsub(Io,panel,de1,de2)
  OPEN(35,FILE=outfile6,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO J=1,Ny  
  DO I=1,Nx  
     WRITE(35, REC=i+(j-1)*Nx) panel(i,j,:)
  END DO
  END DO
  CLOSE(35)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM model3D
