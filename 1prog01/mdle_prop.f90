MODULE mdle_prop

CONTAINS
  SUBROUTINE prop_anisot3D(Vx,Vy,Vz,Sxx,Sxy,Sxz,Syy,Syz,Szz,&
       C11,C12,C13,C33,C44a,C44b,C66,rhox,rhoy,rhoz,dx,dt,Nx,Ny,Nz)
    IMPLICIT NONE
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: Vx,Vz,Sxx,Sxz,Szz
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: Vy,Syy,Sxy,Syz
    REAL, DIMENSION(:,:,:), INTENT(IN)    :: C11,C13,C33,C44a,C44b
    REAL, DIMENSION(:,:,:), INTENT(IN)    :: C12,C66
    REAL, DIMENSION(:,:,:), INTENT(IN)    :: rhox,rhoz
    REAL, DIMENSION(:,:,:), INTENT(IN)    :: rhoy
    REAL,    INTENT(IN) :: dx,dt
    INTEGER, INTENT(IN) :: Nx,Ny,Nz
    REAL :: aux1,aux2,aux3,aux4,aux5,aux6,aux
    INTEGER :: i,j,k
    REAL :: fr1,fr2,ldx

    fr1 = -1./24; fr2=9./8; ldx=1./dx

    ! Space derivatives of stresses => velocities
    DO k=3,Nz-2
       DO j=3,Ny-2
          DO i=3,Nx-2

             aux1 = (fr1)*(Sxx(i+1,j,k)-Sxx(i-2,j,k)) + &
                  (fr2)*(Sxx(i,j,k)-Sxx(i-1,j,k))
             aux1 = aux1*ldx
             aux2 = (fr1)*(Sxy(i,j+1,k)-Sxy(i,j-2,k)) + &
                  (fr2)*(Sxy(i,j,k)-Sxy(i,j-1,k))
             aux2 = aux2*ldx
             aux3 = (fr1)*(Sxz(i,j,k+1)-Sxz(i,j,k-2)) + &
                  (fr2)*(Sxz(i,j,k)-Sxz(i,j,k-1))
             aux3 = aux3*ldx
             aux = dt*rhox(i,j,k)*(aux1+aux2+aux3)
             Vx(i,j,k) = aux + Vx(i,j,k)

             aux1 = (fr1)*(Sxy(i+2,j,k)-Sxy(i-1,j,k)) + &
                  (fr2)*(Sxy(i+1,j,k)-Sxy(i,j,k))
             aux1 = aux1*ldx
             aux2 = (fr1)*(Syy(i,j+2,k)-Syy(i,j-1,k)) + &
                  (fr2)*(Syy(i,j+1,k)-Syy(i,j,k))
             aux2 = aux2*ldx
             aux3 = (fr1)*(Syz(i,j,k+1)-Syz(i,j,k-2)) + &
                  (fr2)*(Syz(i,j,k)-Syz(i,j,k-1))
             aux3 = aux3*ldx
             aux = dt*rhoy(i,j,k)*(aux1+aux2+aux3)
             Vy(i,j,k) = aux + Vy(i,j,k)

             aux1 = (fr1)*(Sxz(i+2,j,k)-Sxz(i-1,j,k)) + &
                  (fr2)*(Sxz(i+1,j,k)-Sxz(i,j,k))
             aux1 = aux1*ldx
             aux2 = (fr1)*(Syz(i,j+1,k)-Syz(i,j-2,k)) + &
                  (fr2)*(Syz(i,j,k)-Syz(i,j-1,k))
             aux2 = aux2*ldx
             aux3 = (fr1)*(Szz(i,j,k+2)-Szz(i,j,k-1)) + &
                  (fr2)*(Szz(i,j,k+1)-Szz(i,j,k))
             aux3 = aux3*ldx
             aux = dt*rhoz(i,j,k)*(aux1+aux2+aux3)
             Vz(i,j,k) = aux + Vz(i,j,k)

          END DO
       END DO
    END DO

    ! Space derivatives of velocities => stresses
    DO k=3,Nz-2
       DO j=3,Ny-2
          DO i=3,Nx-2

             aux1 = (fr1)*(Vx(i+2,j,k)-Vx(i-1,j,k)) + &
                  (fr2)*(Vx(i+1,j,k)-Vx(i,j,k)); aux4=aux1
             aux1 = C11(i,j,k)*aux1*ldx
             aux2 = (fr1)*(Vy(i,j+1,k)-Vy(i,j-2,k)) + &
                  (fr2)*(Vy(i,j,k)-Vy(i,j-1,k)); aux5=aux2
             aux2 = C12(i,j,k)*aux2*ldx
             aux3 = (fr1)*(Vz(i,j,k+1)-Vz(i,j,k-2)) + &
                  (fr2)*(Vz(i,j,k)-Vz(i,j,k-1)); aux6=aux3
             aux3 = C13(i,j,k)*aux3*ldx
             aux = dt*(aux1+aux2+aux3)
             Sxx(i,j,k) = aux + Sxx(i,j,k)

             aux1 = C12(i,j,k)*aux4*ldx
             aux2 = C11(i,j,k)*aux5*ldx
             aux3 = C13(i,j,k)*aux6*ldx
             aux = dt*(aux1+aux2+aux3)
             Syy(i,j,k) = aux + Syy(i,j,k)

             aux1 = C13(i,j,k)*aux4*ldx
             aux2 = C13(i,j,k)*aux5*ldx
             aux3 = C33(i,j,k)*aux6*ldx
             aux = dt*(aux1+aux2+aux3)
             Szz(i,j,k) = aux + Szz(i,j,k)

             aux1 = (fr1)*(Vz(i,j+2,k)-Vz(i,j-1,k)) + &
                  (fr2)*(Vz(i,j+1,k)-Vz(i,j,k))
             aux1 = C44b(i,j,k)*aux1*ldx
             aux2 = (fr1)*(Vy(i,j,k+2)-Vy(i,j,k-1)) + &
                  (fr2)*(Vy(i,j,k+1)-Vy(i,j,k))
             aux2 = C44b(i,j,k)*aux2*ldx
             aux = dt*(aux1+aux2)
             Syz(i,j,k) = aux + Syz(i,j,k)

             aux1 = (fr1)*(Vz(i+1,j,k)-Vz(i-2,j,k)) + &
                  (fr2)*(Vz(i,j,k)-Vz(i-1,j,k))
             aux1 = C44a(i,j,k)*aux1*ldx
             aux2 = (fr1)*(Vx(i,j,k+2)-Vx(i,j,k-1)) + &
                  (fr2)*(Vx(i,j,k+1)-Vx(i,j,k))
             aux2 = C44a(i,j,k)*aux2*ldx
             aux = dt*(aux1+aux2)
             Sxz(i,j,k) = aux + Sxz(i,j,k)

             aux1 = (fr1)*(Vy(i+1,j,k)-Vy(i-2,j,k)) + &
                  (fr2)*(Vy(i,j,k)-Vy(i-1,j,k))
             aux1 = C66(i,j,k)*aux1*ldx
             aux2 = (fr1)*(Vx(i,j+2,k)-Vx(i,j-1,k)) + &
                  (fr2)*(Vx(i,j+1,k)-Vx(i,j,k))
             aux2 = C66(i,j,k)*aux2*ldx
             aux = dt*(aux1+aux2)
             Sxy(i,j,k) = aux + Sxy(i,j,k)

          END DO
       END DO
    END DO


  END SUBROUTINE prop_anisot3D

END MODULE mdle_prop
