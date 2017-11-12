MODULE mdle_io_utils

CONTAINS

  SUBROUTINE inputdata(Nx,Ny,Nz,dx,dt,fpeak,itmax,lx,ly,lz,r,&
       nphones,npmin_x,npmin_y,npmin_z,dnp_x,dnp_y,dnp_z,&
       nsnaps,snapmin,dsnap,iy_snap,nb,F,fsf)

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: Nx,Ny,Nz,itmax,lx,ly,lz
    REAL,    INTENT(OUT) :: dx,dt,fpeak,r
    INTEGER, INTENT(OUT) :: nphones,npmin_x,npmin_y,npmin_z,&
         dnp_x,dnp_y,dnp_z,iy_snap
    INTEGER, INTENT(OUT) :: nsnaps,snapmin,dsnap
    INTEGER, INTENT(OUT) :: nb,fsf
    REAL,    INTENT(OUT) :: F
 
    OPEN(30,FILE="data.dat",STATUS='UNKNOWN',ACTION='READ')
    READ(30,'(t10,i10)') Nx
    READ(30,'(t10,i10)') Ny
    READ(30,'(t10,i10)') Nz
    READ(30,'(t10,f10.4)') dx
    READ(30,'(t10,f10.8)') dt
    READ(30,'(t10,f10.4)') fpeak
    READ(30,'(t10,i10)') itmax
    READ(30,'(t10,i10)') lx
    READ(30,'(t10,i10)') ly
    READ(30,'(t10,i10)') lz
    READ(30,'(t10,f10.4)') r
    READ(30,'(t10,i10)') nb
    READ(30,'(t10,f10.4)') F
    READ(30,'(t10,i10)') fsf
    READ(30,'(t10,i10)') nphones
    READ(30,'(t10,i10)') npmin_x
    READ(30,'(t10,i10)') npmin_y
    READ(30,'(t10,i10)') npmin_z
    READ(30,'(t10,i10)') dnp_x
    READ(30,'(t10,i10)') dnp_y
    READ(30,'(t10,i10)') dnp_z
    READ(30,'(t10,i10)') nsnaps
    READ(30,'(t10,i10)') snapmin
    READ(30,'(t10,i10)') dsnap
    READ(30,'(t10,i10)') iy_snap

  END SUBROUTINE inputdata


  SUBROUTINE inputmodel(Nx,Ny,Nz,C11,C13,C33,C44,C66,rhox)
    IMPLICIT NONE
    REAL, INTENT(OUT), DIMENSION(:,:,:) :: C11,C13,C33,C44,C66,rhox
    INTEGER, INTENT(IN) :: Nx,Ny,Nz
    INTEGER :: i,j

    OPEN(30,FILE="C11.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx
       DO j=1,Ny
          READ(30, REC=i+(j-1)*Nx) C11(i,j,:)
       END DO
    END DO
    CLOSE(30)
    OPEN(30,FILE="C13.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx
       DO j=1,Ny
          READ(30, REC=i+(j-1)*Nx) C13(i,j,:)
       END DO
    END DO
    CLOSE(30)
    OPEN(30,FILE="C33.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx
       DO j=1,Ny
          READ(30, REC=i+(j-1)*Nx) C33(i,j,:)
       END DO
    END DO
    CLOSE(30)
    OPEN(30,FILE="C44.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx
       DO j=1,Ny
          READ(30, REC=i+(j-1)*Nx) C44(i,j,:)
       END DO
    END DO
    CLOSE(30)
    OPEN(30,FILE="C66.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx
       DO j=1,Ny
          READ(30, REC=i+(j-1)*Nx) C66(i,j,:)
       END DO
    END DO
    CLOSE(30)
    OPEN(30,FILE="rhox.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx
       DO j=1,Ny
          READ(30, REC=i+(j-1)*Nx) rhox(i,j,:)
       END DO
    END DO
    CLOSE(30)

  END SUBROUTINE inputmodel



  SUBROUTINE save_shotgather_n_snapshots3D(csg_vx,csg_vy,csg_vz,Vx,Vy,Vz,&
       nphones,npmin_x,npmin_y,npmin_z,dnp_x,dnp_y,dnp_z,&
       isnap,nsnaps,snapmin,dsnap,n1,n2,n3,n4,n5,n6,Nx,Ny,&
       Nz,it,itmax,iy_snap)
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: csg_vx,csg_vy,csg_vz
    REAL, DIMENSION(:,:,:), INTENT(IN)    :: Vx,Vy,Vz
    REAL, DIMENSION(Nx) :: aux1
    REAL, DIMENSION(Nz) :: aux2
    INTEGER, INTENT(IN)    :: nphones
    INTEGER, INTENT(IN)    :: npmin_x,npmin_y,npmin_z,dnp_x,dnp_y,dnp_z
    INTEGER, INTENT(INOUT) :: isnap
    INTEGER, INTENT(IN)    :: nsnaps,snapmin,dsnap
    INTEGER, INTENT(IN)    :: n1,n2,n3,n4,n5,n6
    INTEGER, INTENT(IN)    :: Nx,Ny,Nz,it,itmax,iy_snap
    INTEGER :: i,ix,iy,iz

    do i=1,nphones
       ix = (i-1)*dnp_x + npmin_x
       iy = (i-1)*dnp_y + npmin_y
       iz = (i-1)*dnp_z + npmin_z
       csg_vx(it,i) = Vx(ix,iy,iz)
       csg_vy(it,i) = Vy(ix,iy,iz)
       csg_vz(it,i) = Vz(ix,iy,iz)
    end do

    if (it == itmax) then  
       PRINT*,'Saving CSG'
       do i=1,nphones
          WRITE(n1, REC=i) csg_vx(:,i)
          WRITE(n2, REC=i) csg_vy(:,i)
          WRITE(n3, REC=i) csg_vz(:,i)
       end do
    end if

    if (it == (isnap * dsnap) + snapmin) then
       isnap=isnap+1
       if (isnap <= nsnaps) then
          PRINT*, 'Saving snapshot ',isnap,'/',nsnaps
          do i=1,Nx
             WRITE(n4, REC=((isnap -1)*Nx + i)) Vx(i,iy_snap,:)
             WRITE(n5, REC=((isnap -1)*Nx + i)) Vy(i,iy_snap,:)
             WRITE(n6, REC=((isnap -1)*Nx + i)) Vz(i,iy_snap,:)
          end do
       end if
    end if

  END SUBROUTINE save_shotgather_n_snapshots3D

END MODULE mdle_io_utils


