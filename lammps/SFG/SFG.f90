!Set of subroutines used to analyse the lipid trajectory

MODULE histogram
  IMPLICIT NONE 
  REAL*8          , ALLOCATABLE            :: oh_dist(:,:,:)
  REAL*8          , ALLOCATABLE            :: oh_vel(:,:,:)
  REAL*8          , ALLOCATABLE            :: autocorr(:,:)
  REAL*8          , ALLOCATABLE            :: layerspf(:)
END MODULE histogram

PROGRAM SFG
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_CUSTOM_POS_VEL
    CALL SET_CENTER_OF_MASS_TO_ZERO
    CALL OH_DISTRIBUTION (frame)
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
                         nequil,  nwater, dt, layers, &
                         atype, moltype
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
  ALLOCATE(layers(nlayers+1))
  READ(1,*) layers
  CLOSE(1)

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(vel(natoms,3)) ! velocities 
  ALLOCATE(coarse(nframes,2*nwater)); coarse=0
  ALLOCATE(oh_dist(3,2*nwater,nframes)); oh_dist = 0.d0
  ALLOCATE(oh_vel(3,2*nwater,nframes)); oh_vel = 0.d0
  ALLOCATE(autocorr(nlayers,nframes)); autocorr = 0
  ALLOCATE(moltype(natoms))
  ALLOCATE(atype(natoms))

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = vel_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE OH_DISTRIBUTION (frame)
  !Generate the OH distribution for each water in a specific frame
  USE histogram,  ONLY : oh_dist, oh_vel
  USE parameters, ONLY : natoms, vel, coarse, moltype
  IMPLICIT NONE
  INTEGER                    :: frame, i, ih
  INTEGER                    :: ind_OH
  INTEGER                    :: ind_layer
  INTEGER                    :: ind_H(2)
  LOGICAL                    :: is_oxygen
  REAL*8                     :: oh(3,2)
  LOGICAL                    :: is_water_oxygen

  DO i = 1,natoms

    !Calculate the OH vector:
    IF ( is_water_oxygen(i) ) THEN

      ind_H = 0
      ind_OH = 2*(moltype(i)-2) !Assuming water indexes start from 2
      CALL OH_VECT( i, oh, ind_H )
   
      !Max two OH vectors per O atom, taking PBC into account 
      DO ih = 1,2

        ind_OH = ind_OH + 1
        oh_dist (:,ind_OH,frame) = oh(:,ih)  
        oh_vel (:,ind_OH,frame) = vel(ind_H(ih),:) - vel(i,:) 
        coarse(frame,ind_OH) = ind_layer(i)

      END DO

    END IF

  ENDDO
 
END SUBROUTINE OH_DISTRIBUTION

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
        IF (layer.eq.4) CYCLE
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

SUBROUTINE OH_VECT(O, OH, ind_H)
  ! Normalized vector from atom O to atom H
  USE parameters, ONLY : natoms, box, pos, moltype
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: O
  INTEGER                              :: ind_H(2)
  REAL*8          , INTENT(INOUT)      :: OH(3,2)
  REAL*8                               :: ohtmp(3)
  INTEGER                              :: iat, ipol, ih

  ih=0
 
  DO iat=1,natoms
 
    IF ( ( iat.ne.O ) .and. ( moltype(iat).eq.moltype(O) ) ) THEN

      DO ipol=1,3
    
        ohtmp(ipol) = pos(iat,ipol) - pos(O,ipol)
        ohtmp(ipol) = ohtmp(ipol) - nint( ohtmp(ipol)/box(ipol) ) * box(ipol)
    
      END DO
    
      ih = ih + 1
      OH(:,ih) = ohtmp / norm2(ohtmp)
      ind_H(ih) = iat
 
    END IF

  END DO

  IF ( ih.ne.2 ) THEN
    PRINT *, 'Wrong number of OH bonds per water. STOP!!!', ih
    STOP
  END IF

  IF ( ind_H(1) > ind_H(2) ) THEN
    ohtmp=OH(:,2)
    OH(:,2)=OH(:,1)
    OH(:,1)=ohtmp
    iat=ind_H(1)
    ind_H(1)=ind_H(2)
    ind_H(2)=iat
  END IF

END SUBROUTINE OH_VECT

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
