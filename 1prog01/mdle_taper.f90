MODULE mdle_taper

CONTAINS

  SUBROUTINE bt_exp_create(taper,nb,F)
    IMPLICIT NONE
    REAL :: F
    INTEGER :: nb
    REAL, DIMENSION(nb) :: taper
    INTEGER :: i

    DO i=1,nb
       taper(i) = exp( -(F*(REAL(nb) - REAL(i) ))**2 )
    END DO
  END SUBROUTINE bt_exp_create

  SUBROUTINE bt_apply_multiple(Vx,Vy,Vz,Sxx,Sxy,Sxz,&
       Syy,Syz,Szz,Nx,Ny,Nz,nb,taper)
    IMPLICIT NONE
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: Vx,Vy,Vz,Sxx,Sxy,Sxz
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: Syy,Syz,Szz
    INTEGER, INTENT(IN) :: Nx,Ny,Nz,nb
    REAL, DIMENSION(nb), INTENT(IN) :: taper
    INTEGER :: i

    DO i=1,Ny
       CALL bt_apply(Vx(:,i,:),Nx,Nz,nb,taper)
       CALL bt_apply(Vy(:,i,:),Nx,Nz,nb,taper)
       CALL bt_apply(Vz(:,i,:),Nx,Nz,nb,taper)
       CALL bt_apply(Sxx(:,i,:),Nx,Nz,nb,taper)
       CALL bt_apply(Sxy(:,i,:),Nx,Nz,nb,taper)
       CALL bt_apply(Sxz(:,i,:),Nx,Nz,nb,taper)
       CALL bt_apply(Syy(:,i,:),Nx,Nz,nb,taper)
       CALL bt_apply(Syz(:,i,:),Nx,Nz,nb,taper)
       CALL bt_apply(Szz(:,i,:),Nx,Nz,nb,taper)
    END DO
    DO i=1,Nx
       CALL bt_apply(Vx(i,:,:),Nx,Nz,nb,taper)
       CALL bt_apply(Vy(i,:,:),Nx,Nz,nb,taper)
       CALL bt_apply(Vz(i,:,:),Nx,Nz,nb,taper)
       CALL bt_apply(Sxx(i,:,:),Nx,Nz,nb,taper)
       CALL bt_apply(Sxy(i,:,:),Nx,Nz,nb,taper)
       CALL bt_apply(Sxz(i,:,:),Nx,Nz,nb,taper)
       CALL bt_apply(Syy(i,:,:),Nx,Nz,nb,taper)
       CALL bt_apply(Syz(i,:,:),Nx,Nz,nb,taper)
       CALL bt_apply(Szz(i,:,:),Nx,Nz,nb,taper)
    END DO
    DO i=1,Nz
       CALL bt_apply(Vx(:,:,i),Nx,Nz,nb,taper)
       CALL bt_apply(Vy(:,:,i),Nx,Nz,nb,taper)
       CALL bt_apply(Vz(:,:,i),Nx,Nz,nb,taper)
       CALL bt_apply(Sxx(:,:,i),Nx,Nz,nb,taper)
       CALL bt_apply(Sxy(:,:,i),Nx,Nz,nb,taper)
       CALL bt_apply(Sxz(:,:,i),Nx,Nz,nb,taper)
       CALL bt_apply(Syy(:,:,i),Nx,Nz,nb,taper)
       CALL bt_apply(Syz(:,:,i),Nx,Nz,nb,taper)
       CALL bt_apply(Szz(:,:,i),Nx,Nz,nb,taper)
    END DO

  END SUBROUTINE bt_apply_multiple


  SUBROUTINE bt_apply(pp,Nx,Nz,nb,taper)
    IMPLICIT NONE
    INTEGER :: NX,NZ,nb
    REAL, DIMENSION(nb) :: taper
    REAL, DIMENSION(NX,NZ) :: pp
    INTEGER :: i

    DO i=1,NZ
       pp(1:nb,i) = pp(1:nb,i) * taper
       pp(NX:NX-nb+1:-1,i) =  pp(NX:NX-nb+1:-1,i) * taper
    END DO

    DO i=1,NX
       pp(i,1:nb:1) = pp(i,1:nb:1) * taper
       pp(i,NZ:NZ-nb+1:-1) =  pp(i,NZ:NZ-nb+1:-1) * taper  
    END DO
  END SUBROUTINE bt_apply

END MODULE mdle_taper



