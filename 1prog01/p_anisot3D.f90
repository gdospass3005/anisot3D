PROGRAM anisot

  ! Modeling seismic waves in 3-D TI media

  USE mdle_source
  USE mdle_prop
  USE mdle_io_utils                                                  
  USE mdle_taper

  IMPLICIT NONE
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: Vx,Vz,Sxx,Sxz,Szz
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: Vy,Syy,Sxy,Syz

  INTEGER :: Nx,Ny,Nz
  REAL    :: dx,dt,fpeak
  INTEGER :: itmax
  INTEGER :: lx, lz
  INTEGER :: ly
  REAL    :: r

  REAL, DIMENSION(:,:,:), ALLOCATABLE :: C11,C13,C33,C44,C44a,C44b,C12,C66
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: rhox,rhoy,rhoz

  !REAL, DIMENSION(:,:,:), ALLOCATABLE :: source
  REAL, DIMENSION(:),   ALLOCATABLE :: trs
  REAL    :: tdelay

  REAL :: F
  REAL, DIMENSION(:),   ALLOCATABLE :: taper
  INTEGER :: nb,fsf

  REAL    :: t
  INTEGER :: it,i,j,k

  INTEGER :: iy_csg
  REAL, DIMENSION(:,:), ALLOCATABLE :: csg_vx,csg_vz
  REAL, DIMENSION(:,:), ALLOCATABLE :: csg_vy
  INTEGER :: nphones,npmin_x,npmin_z,dnp_x,dnp_z
  INTEGER :: npmin_y,dnp_y
  INTEGER :: n1,n2,n3,n4,n5,n6
  INTEGER :: isnap,nsnaps,snapmin,dsnap,iy_snap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PRINT*,'anisot3D - Modeling seismic waves in 3-D anisotropic media'

  CALL inputdata(Nx,Ny,Nz,dx,dt,fpeak,itmax,lx,ly,lz,r,&
       nphones,npmin_x,npmin_y,npmin_z,dnp_x,dnp_y,dnp_z,&
       nsnaps,snapmin,dsnap,iy_snap,nb,F,fsf)

  ALLOCATE(taper(nb))
  CALL bt_exp_create(taper,nb,F)

  ALLOCATE(C11(Nx,Ny,Nz),C13(Nx,Ny,Nz),C33(Nx,Ny,Nz),C44(Nx,Ny,Nz))
  ALLOCATE(C12(Nx,Ny,Nz),C66(Nx,Ny,Nz))
  ALLOCATE(rhox(Nx,Ny,Nz),rhoz(Nx,Ny,Nz))
  ALLOCATE(rhoy(Nx,Ny,Nz))

  CALL inputmodel(Nx,Ny,Nz,C11,C13,C33,C44,C66,rhox)

  C11=C11*1e10;   C33=C33*1e10;   C44=C44*1e10;   C13=C13*1e10; 
  C66=C66*1e10; 
  C12 = C11 - 2.0*C66

  PRINT*,'MAXVAL(C11)',MAXVAL(C11)
  PRINT*,'MAXVAL(C33)',MAXVAL(C33)
  PRINT*,'MAXVAL(C44)',MAXVAL(C44)
  PRINT*,'MAXVAL(C12)',MAXVAL(C12)
  PRINT*,'MAXVAL(C13)',MAXVAL(C13)


  do k=1,Nz
     do j=1,Ny
        do i=1,Nx-1
           C11(i,j,k) = 0.5*(C11(i,j,k) + C11(i+1,j,k))
           C12(i,j,k) = 0.5*(C12(i,j,k) + C12(i+1,j,k))
           C13(i,j,k) = 0.5*(C13(i,j,k) + C13(i+1,j,k))
           C33(i,j,k) = 0.5*(C33(i,j,k) + C33(i+1,j,k))
        end do
     end do
  end do
  do k=1,Nz
     do j=1,Ny-1
        do i=1,Nx
           C66(i,j,k) = 0.5*(C66(i,j,k) + C66(i,j+1,k))
        end do
     end do
  end do
  ALLOCATE(C44a(Nx,Ny,Nz),C44b(Nx,Ny,Nz))
  C44a=C44; C44b=C44
  do k=1,Nz
     do j=1,Ny-1
        do i=1,Nx
           C44a(i,j,k) = 0.5*(C44(i,j,k) + C44(i,j+1,k))
        end do
     end do
  end do
  do k=1,Nz-1
     do j=1,Ny-1
        do i=1,Nx-1
           C44b(i,j,k) = (1./8)*(C44(i,j,k) + C44(i,j,k+1) &
                + C44(i,j+1,k) + C44(i,j+1,k+1) &
                + C44(i+1,j,k) + C44(i+1,j,k+1) &
                + C44(i+1,j+1,k) + C44(i+1,j+1,k+1) )
        end do
     end do
  end do
  DEALLOCATE(C44)
  rhoy=rhox
  do k=1,Nz
     do j=1,Ny-1
        do i=1,Nx-1
           rhoy(i,j,k) = 0.25*(rhox(i,j,k) + rhox(i+1,j,k) &
                + rhox(i,j+1,k) + rhox(i+1,j+1,k))
        end do
     end do
  end do
  rhoz=rhox
  do k=1,Nz-1
     do j=1,Ny
        do i=1,Nx-1
           rhoz(i,j,k) = 0.25*(rhox(i,j,k) + rhox(i+1,j,k) &
                + rhox(i,j,k+1) + rhox(i+1,j,k+1))
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(trs(itmax))
  CALL source_init(trs,dt,itmax,fpeak,tdelay)
  OPEN(76,FILE="source.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  WRITE(76, REC=1) trs(:)

  PRINT*,'tdelay',tdelay
  PRINT*,'amostras da fonte:',FLOOR(2*tdelay/dt)

  ALLOCATE(csg_vx(itmax,nphones),csg_vz(itmax,nphones))
  ALLOCATE(csg_vy(itmax,nphones))

  n1=31; n2=32; n3=33; n4=34; n5=35; n6=36
  OPEN(n1,FILE="csg_Vx.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  OPEN(n2,FILE="csg_Vy.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  OPEN(n3,FILE="csg_Vz.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  OPEN(n4,FILE="snapshots_Vx.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  OPEN(n5,FILE="snapshots_Vy.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  OPEN(n6,FILE="snapshots_Vz.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(Vx(Nx,Ny,Nz),Vz(Nx,Ny,Nz),Sxx(Nx,Ny,Nz),Sxz(Nx,Ny,Nz),&
       Szz(Nx,Ny,Nz))
  ALLOCATE(Vy(Nx,Ny,Nz),Syy(Nx,Ny,Nz),Sxy(Nx,Ny,Nz),Syz(Nx,Ny,Nz))

  Sxx=0.; Sxz=0.; Szz=0.; Vx=0.; Vz=0.
  Syy=0.; Sxy=0.; Syz=0.; Vy=0.
  rhox=1./rhox; rhoz=1./rhoz
  rhoy=1./rhoy
  isnap=0

  ! Beggining time loop
  DO it=1,itmax
     t=it*dt
     OPEN(77,FILE="status.txt",STATUS='UNKNOWN',ACTION='WRITE')
     WRITE(77,*) 'anisot - Modeling seismic waves in 3-D anisotropic media'
     WRITE(77,*) 'Iteracao',it,'/',itmax,' Vx(lx+5,ly+5,lz+5)=',&
          Vx(lx+5,ly+5,lz+5)
     CLOSE(77)

     ! Insert source function 
     Sxx(lx,ly,lz) = Sxx(lx,ly,lz) + trs(it)
     Syy(lx,ly,lz) = Syy(lx,ly,lz) + trs(it)
     Szz(lx,ly,lz) = Szz(lx,ly,lz) + trs(it)

     ! Do one time step
     CALL prop_anisot3D(Vx,Vy,Vz,Sxx,Sxy,Sxz,Syy,Syz,Szz,&
          C11,C12,C13,C33,C44a,C44b,C66,rhox,rhoy,rhoz,dx,dt,Nx,Ny,Nz)

     PRINT*,'it',it,'/',itmax,' Vx(lx+5,ly+5,lz+5)=',Vx(lx+5,ly+5,lz+5)

     ! Output
     CALL save_shotgather_n_snapshots3D(csg_vx,csg_vy,csg_vz,Vx,Vy,Vz,&
          nphones,npmin_x,npmin_y,npmin_z,dnp_x,dnp_y,dnp_z,&
          isnap,nsnaps,snapmin,dsnap,n1,n2,n3,n4,n5,n6,Nx,Ny,&
          Nz,it,itmax,iy_snap)

     ! Boundary taper
     CALL bt_apply_multiple(Vx,Vy,Vz,Sxx,Sxy,Sxz,Syy,Syz,Szz,&
          Nx,Ny,Nz,nb,taper)

  END DO
  PRINT*,'Successful run.'

END PROGRAM anisot
