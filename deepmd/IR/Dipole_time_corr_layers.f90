
MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: dipt(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
END MODULE histogram

PROGRAM VAC 
  USE parameters, ONLY : nframes, coarse
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_RAW_WANNIER 
    CALL READ_RAW_POS_BOX
    !CALL SET_CENTER_OF_MASS_TO_ZERO
    CALL ASSIGN_DIPOLE (frame) 
    CALL COARSE_GRAIN_POS (frame)
    CALL MAKE_NEIGHBOR_LIST (frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL COARSE_GRAIN_HYDROGEN
  CALL DIPOLE_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM VAC

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : dipt, autocorr 
  USE parameters, ONLY : pos, wannier, natoms, nwater, &
                         nframes, coarse, nlayers,&
                         layers, dt
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil
  READ(1,*) nwater, nlayers, dt
  ALLOCATE(layers(nlayers+1))
  READ(1,*) layers
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(wannier(3,4,nwater)) 
  ALLOCATE(dipt(3,nwater,nframes)) ! dipole moment of all water
  ALLOCATE(coarse(nframes,nwater)); coarse = 1
  ALLOCATE(autocorr(nframes,nlayers)); autocorr = 0

  OPEN(unit = 1,file = 'coord.raw')
  OPEN(unit = 2,file = 'box.raw')
  OPEN(unit = 3,file = 'type.raw')
  OPEN(unit = 4,file = 'wannier.raw')
  
END SUBROUTINE INITIALIZE

SUBROUTINE ASSIGN_DIPOLE(frame)
  USE histogram, ONLY : dipt
  USE parameters, ONLY : wannier, natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, iw, iwfc
  REAL*8                               :: d(4), posOw(3)
  
  iw = 0
  dipt(:,:,frame) = 0

  DO iat=1,natoms

    !Electronic part
    IF ( is_water_oxygen(iat) ) THEN
  
      iw = iw + 1
      DO iwfc = 1, 4
        !Assuming Wannier coordinates are relative to oxygen atom
        dipt(:,iw,frame) = dipt(:,iw,frame) - 2*wannier(:,iwfc,iw)
      END DO
  
      posOw = pos(:,iat)

    ELSE

      CALL DISTANCE_VECTOR( posOw, pos(:,iat), d ) 
      !Ionic part
      dipt(:,iw,frame) = dipt(:,iw,frame) + d(1:3) 

    END IF

  END DO

END SUBROUTINE ASSIGN_VEL

SUBROUTINE DIPOLE_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : dipt, autocorr
  USE parameters, ONLY : nframes, natoms, iwater, coarse, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_dip(nlayers)
  DOUBLE PRECISION           :: dotprod, surv_prob(nlayers)
  INTEGER                    :: frame1, frame2, i_at, i_coord
  INTEGER                    :: layer, i_bin
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers) 
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,nframes-1

    surv_prob = 0.d0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,20

      counter_in = 0
      counter_all = 0
      sum_dip = 0.D0

!$omp parallel do private(layer,dotprod) reduction(+:counter_in, counter_all, sum_dip)
      DO i_at = 1,nwater

        layer = coarse(frame2,iwater(i_at,frame2))
        counter_all(layer) = counter_all(layer) + 1

        IF ( in_layer(i_at,frame2,frame2+frame1) ) THEN

          counter_in(layer) = counter_in(layer) + 1

          dotprod = sum( dipt(:,i_at,frame2) * &
                         dipt(:,i_at, frame2+frame1) )

          sum_dip(layer) = sum_dip(layer) + dotprod

        ENDIF
      ENDDO
!$omp end parallel do

      DO i_bin = 1,nlayers
        IF ( counter_all(i_bin) /= 0 ) THEN
          surv_prob(i_bin) = surv_prob(i_bin) &
                           + dble(counter_in(i_bin)) / dble(counter_all(i_bin))
          !Averaged correlation function
          autocorr(frame1+1,i_bin) = autocorr(frame1+1,i_bin) + &
                                sum_dip(i_bin) / dble( counter_all(i_bin) ) 
        END IF
      END DO

    ENDDO

    surv_prob = surv_prob / dble(frame2-1) 

    DO i_bin = 1, nlayers

      IF ( surv_prob(i_bin) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(frame1+1,i_bin) = autocorr(frame1+1,i_bin) / dble( frame2-1 ) / surv_prob(i_bin) 
      END IF

    END DO

  ENDDO

END SUBROUTINE DIPOLE_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,   ONLY : autocorr
  USE parameters, ONLY : box, dt, nlayers, nframes
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION,PARAMETER :: au_to_fs = 2.418884326509D-2

  OPEN(unit = 3,file = 'vac.dat')

  DO i = 0,nframes-1

    write(3,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*dt*au_to_fs, autocorr(i+1,:) / autocorr(1,:) 

  ENDDO

  CLOSE(3)

END SUBROUTINE PRINT_RESULTS

