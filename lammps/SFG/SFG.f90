!Set of subroutines used to analyse the lipid trajectory

MODULE histogram
  IMPLICIT NONE 
  REAL*8          , ALLOCATABLE            :: oh_dist(:,:,:)
  REAL*8          , ALLOCATABLE            :: oh_vel(:,:,:)
  REAL*8          , ALLOCATABLE            :: autocorr(:,:)
END MODULE histogram

PROGRAM SFG
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  CALL REMOVE_EQUIL_VEL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL READ_VELOCITIES
    CALL COARSE_GRAIN_POS(frame)
    CALL OH_DISTRIBUTION(frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL SSVAC 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM SFG

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : autocorr, oh_dist, oh_vel
  USE parameters, ONLY : natoms, nframes, nwater, &
                         pos, vel, coarse, nlayers, &
                         nequil, dt, layers, &
                         atypeO, atypeH, atype
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, vel_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, vel_file)
  CALL getarg(3, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil
  READ(1,*) nwater, nlayers, dt
  READ(1,*) atypeO, atypeH
  ALLOCATE(layers(nlayers+1))
  READ(1,*) layers
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(vel(3,natoms)) ! velocities 
  ALLOCATE(coarse(nframes,2*nwater)); coarse=1
  ALLOCATE(oh_dist(3,2*nwater,nframes)); oh_dist = 0.d0
  ALLOCATE(oh_vel(3,2*nwater,nframes)); oh_vel = 0.d0
  ALLOCATE(autocorr(nlayers,nframes)); autocorr = 0
  ALLOCATE(atype(natoms))

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = vel_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE SSVAC 
  !Surface specific velocity-velocity correlation function
  !Check eq. 9 of j. chem. phys 143, 124702 (2015)
  USE histogram,  ONLY : oh_dist, oh_vel, autocorr
  USE parameters, ONLY : nframes, nwater, coarse, nlayers
  IMPLICIT NONE
  REAL*8                     :: sum_dip(nlayers), surv_prob(nlayers)
  REAL*8                     :: second_leg, dotprod 
  INTEGER                    :: frame1, frame2, iw, ipol
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers), ih 
  INTEGER                    :: layer
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,(nframes-1)/5.

    surv_prob = 0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,100

      counter_in = 0
      counter_all = 0
      sum_dip = 0.D0

!$omp parallel do private(layer,dotprod) reduction(+:counter_in, sum_dip, counter_all)
      DO iw = 1,2*nwater

        layer = coarse(frame2,iw)
        counter_all(layer) = counter_all(layer) + 1

        !If the particle dont leave the window during the interval frame1
        IF ( in_layer(iw,frame2,frame2+frame1,100) ) THEN

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
                                sum_dip(layer)  
        END IF
      END DO

    ENDDO

    DO layer = 1, nlayers
      IF ( surv_prob(layer) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(layer,frame1+1) = autocorr(layer,frame1+1) / dble( frame2-1 ) * 100
      END IF
    END DO

  ENDDO

END SUBROUTINE SSVAC 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : autocorr
  USE parameters, ONLY : nlayers, nframes, dt
  IMPLICIT NONE
  INTEGER                    :: i, j

  OPEN( unit = 2, file = "ssvac.dat" )

  DO i = 0,nframes-1

    write(2,fmt = '(E13.6,3X,*(E14.7,3X))'), dble(i)*dt, &
         ( autocorr(j,i+1) , j=1,nlayers) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE OH_DISTRIBUTION (frame)
  !Compute OH_distance and OH velocity for each H atom
  USE histogram,  ONLY : oh_dist, oh_vel
  USE parameters, ONLY : natoms, vel, atypeO, layers
  IMPLICIT NONE
  INTEGER                    :: frame, i, ih
  INTEGER                    :: ind_OH, indO
  REAL*8                     :: d_oh(3)
  LOGICAL                    :: is_hydrogen

  ind_OH=0
  DO i = 1,natoms
    IF ( is_hydrogen(i) ) THEN

      ind_OH = ind_OH + 1
      CALL MINDISTANCE_VECTOR_IND( i, atypeO, d_oh, indO )
      oh_dist (:,ind_OH,frame) = d_oh / norm2(d_oh)
      oh_vel (:,ind_OH,frame) = vel(:,i) - vel(:,indO)

    END IF

  ENDDO
 
END SUBROUTINE OH_DISTRIBUTION

INTEGER FUNCTION ind_layer( ind_O )
  !Layer index for prefactor. ind_O is the index of water oxygen in pos
  USE parameters, ONLY : pos, box, layers, nlayers
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind_O
  INTEGER                    :: il
  REAL*8                     :: z

  z = pos(ind_O,3)
  z = z - nint(z/box(3))*box(3)
  z = z + box(3)/2.d0

  DO il = 1,nlayers

    if ( z <= layers(il+1) ) THEN
      ind_layer = il
      RETURN
    endif

  END DO

  ind_layer = 0

END FUNCTION ind_layer
