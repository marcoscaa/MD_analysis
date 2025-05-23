!Set of subroutines used to analyse the lipid trajectory

MODULE parameters
  IMPLICIT NONE 
  REAL*8, ALLOCATABLE                      :: pos(:,:,:), oh_vel(:,:,:)
  REAL*8, ALLOCATABLE                      :: vel(:,:,:), oh_dist(:,:,:) 
  REAL*8, ALLOCATABLE                      :: autocorr(:,:)
  REAL*8                                   :: box(3)
  INTEGER*8, ALLOCATABLE                   :: coarse(:,:)
  INTEGER,PARAMETER                        :: nframes = 20000, nlayers = 2
  REAL*8,PARAMETER                         :: dt = 0.5 ! in fs
  INTEGER                                  :: natoms, nwater

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  DO frame = 1,nframes
    CALL READ_FRAME (frame)
    CALL COARSE_GRAIN (frame)
    CALL COMPUTE_OH_VELOCITIES (frame )
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL SSVAC 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file
  CHARACTER(5)               :: moltype,atomtype
  INTEGER                    :: junk1

  ! init
  nwater = 0

  !User should provide filename and index_file as arguments 
  CALL getarg(1, file_name)

  !Assuming constant number of particles during the simulation
  !Getting the number of atoms only once
  OPEN(unit = 1, file = file_name)
  READ(1,*)
  READ(1,*), natoms
  nwater=natoms/3
  CLOSE(1)

  !!!!!!!Allocating all the arrays to be used
  ALLOCATE(pos(3,3,nwater)) ! coordinates, atomtype (O, H1 and H2), n_water
  ALLOCATE(vel(3,3,nwater)) ! coordinates, atomtype (O, H1 and H2), n_water
  ALLOCATE(coarse(nframes,2*nwater)); coarse = 0
  ALLOCATE(oh_vel(3,2*nwater,nframes)); oh_vel = 0.d0
  ALLOCATE(oh_dist(3,2*nwater,nframes)); oh_dist = 0.d0
  ALLOCATE(autocorr(nlayers,nframes)); autocorr = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  OPEN(unit = 1,file = file_name)
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME (frame)
  !Subroutine reads gro files
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: junk1, frame, io 
  REAL*8                     :: z_pbc, z
  CHARACTER(5)               :: moltype,atomtype

  !Discarding first two lines
  READ(1,*)
  READ(1,*)

  DO io = 1,nwater
    READ(1, fmt="(i5,2a5,i5,3f8.3,3f8.4)"), junk1, moltype, atomtype, junk1, pos(:,1,io), vel(:,1,io)
    READ(1, fmt="(i5,2a5,i5,3f8.3,3f8.4)"), junk1, moltype, atomtype, junk1, pos(:,2,io), vel(:,2,io)
    READ(1, fmt="(i5,2a5,i5,3f8.3,3f8.4)"), junk1, moltype, atomtype, junk1, pos(:,3,io), vel(:,3,io)
  END DO

  READ(1,*), box(1), box(2), box(3)

END SUBROUTINE

SUBROUTINE COARSE_GRAIN (frame)
  !Create a coarse grain function to determine (with time) if an atom is
  !or not inside a specific window
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: j, frame
  REAL*8                     :: z, z_pbc
  REAL*8                     :: lower, upper 

  upper = -0.5
  lower = -1.5

  DO j = 1, nwater

    !Considering pbc in the z direction
    z = z_pbc( pos(3,1,j), box(3) )

    IF ( z <= lower )  THEN
      coarse(frame,2*j) = 1
      coarse(frame,2*j-1) = 1
    ELSEIF (( z >= upper ) .and. ( z <= 0.3 )) THEN
      coarse(frame,2*j) = 2
      coarse(frame,2*j-1) = 2
    ENDIF

  ENDDO

END SUBROUTINE COARSE_GRAIN 

SUBROUTINE COMPUTE_OH_VELOCITIES (frame)
  USE parameters
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: frame
  INTEGER                :: iw,ih, ipol
  REAL*8                 :: oh_tmp(3)

  do iw=1,nwater

    do ih=1,2

      oh_vel(:,2*iw-ih+1,frame) = vel(:,ih+1,iw) - vel(:,1,iw)

      do ipol = 1,3

        oh_tmp(ipol) = pos(ipol,ih+1,iw) - pos(ipol,1,iw)
        oh_tmp(ipol) = oh_tmp(ipol) - nint(oh_tmp(ipol)/box(ipol))*box(ipol)

      end do

      oh_dist(:,2*iw-ih+1,frame) = oh_tmp/norm2(oh_tmp)

    end do

  end do

END SUBROUTINE COMPUTE_OH_VELOCITIES

SUBROUTINE SSVAC 
  !Surface specific velocity-velocity correlation function
  !Check eq. 9 of j. chem. phys 143, 124702 (2015)
  USE parameters,  ONLY : oh_dist, oh_vel, autocorr, &
                          nframes, nwater, coarse, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_dip(nlayers), surv_prob(nlayers)
  DOUBLE PRECISION           :: dotprod 
  INTEGER                    :: frame1, frame2, iw, ipol
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers), ih 
  INTEGER                    :: layer
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,(nframes-1)/5.

    surv_prob = 0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,10

      counter_in = 0
      counter_all = 0
      sum_dip = 0.D0

!$omp parallel do private(layer,dotprod) reduction(+:counter_in, sum_dip, counter_all)
      DO iw = 1,2*nwater

        layer = coarse(frame2,iw)
        counter_all(layer) = counter_all(layer) + 1

        IF ( layer .eq. 0 ) CYCLE

        !If the particle dont leave the window during the interval frame1
        IF ( in_layer(iw,frame2,frame2+frame1,10) ) THEN

          counter_in(layer) = counter_in(layer) + 1

          dotprod = 0.d0

          DO ipol = 1,3

            dotprod = dotprod + &
                    & oh_vel(ipol, iw, frame2+frame1) * &
                    & oh_dist(ipol, iw, frame2+frame1)
          END DO

          dotprod = dotprod * oh_vel(3, iw, frame2)
          sum_dip(layer) = sum_dip(layer) + dotprod 

        ENDIF

      ENDDO
!$omp end parallel do

      DO layer = 1,nlayers
        IF ( counter_all(layer) /= 0 ) THEN
          surv_prob(layer) = surv_prob(layer) &
                           + dble(counter_in(layer)) / dble(counter_all(layer))
          !Averaged correlation function
          autocorr(layer,frame1+1) = autocorr(layer,frame1+1) + &
                                sum_dip(layer) / dble( counter_all(layer) ) 
        END IF
      END DO

    ENDDO

    surv_prob = 10 * surv_prob / dble(frame2-1)  

    DO layer = 1, nlayers
      IF ( surv_prob(layer) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(layer,frame1+1) = 10 * autocorr(layer,frame1+1) / dble( frame2-1 ) / surv_prob(layer) 
      END IF
    END DO

  ENDDO

END SUBROUTINE SSVAC 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE parameters, ONLY : nlayers, nframes, dt, autocorr
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION           :: n_const(nlayers)

  OPEN( unit = 2, file = "ssvac.dat" )

  DO i = 0,nframes-1

    write(2,fmt = '(E13.6,3X,*(E14.7,3X))'), dble(i)*dt, &
         ( autocorr(j,i+1) / autocorr(j,1) , j=1,nlayers) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

REAL*8 FUNCTION Dist(a, b)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, INTENT(IN)       :: a(3),b(3)
  REAL*8                   :: xyz(3)
  INTEGER                  :: i

  DO i = 1,3

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

REAL*8 FUNCTION z_pbc(z, box_z)
  !Apply PBC for the z coordinate olny
  IMPLICIT NONE
  REAL*8, INTENT(IN)               :: z, box_z

  z_pbc = z - nint(z/box_z)*box_z

END FUNCTION z_pbc

LOGICAL FUNCTION in_layer( iat, ti, tf, sep )
  !Test if atom is inside a layer from time ti to tf
  USE parameters, ONLY : coarse, nlayers
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: iat, ti, tf
  INTEGER                    :: layer, ts
  INTEGER, OPTIONAL          :: sep

  ts = 1
  IF (present(sep)) ts=sep

  IF ( nlayers > 1 ) THEN

    layer = coarse(ti,iat)
    in_layer = ALL( coarse(ti:tf:ts,iat) == layer )

  ELSE 

    in_layer = .true.

  END IF

END FUNCTION in_layer

SUBROUTINE SMOOTH_COARSE_GRAIN
  !Avoid fast jump between layers in the coarse grain function
  USE parameters, ONLY : coarse, nwater, nframes, nlayers
  IMPLICIT NONE
  INTEGER                    :: iw, frame, layer, tout, np
  INTEGER, PARAMETER         :: forgiv=20
  
  IF (nlayers < 2) RETURN

  np = size(coarse,2)

  DO iw = 1, np
  
    frame = 1
    tout  = 0
    layer = coarse(frame,iw)

    DO WHILE ( frame <= nframes ) 

      IF ( coarse(frame,iw) /= layer ) THEN
        tout = tout + 1
      ELSEIF ( tout > 0 ) THEN 
        coarse(frame-tout:frame,iw) = layer
        tout = 0
      END IF
  
      !If the particle really changed layer
      IF ( tout > forgiv ) THEN

        layer = coarse(frame,iw)
        tout = 0

      END IF

      frame = frame + 1

    END DO

  END DO

END SUBROUTINE SMOOTH_COARSE_GRAIN
