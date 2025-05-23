MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: lifetime(:)
  DOUBLE PRECISION                         :: timestep
  INTEGER, ALLOCATABLE                     :: coarse(:,:)
  INTEGER                                  :: ndata, nsteps
END MODULE histogram

PROGRAM OHcorr
  IMPLICIT NONE

  CALL INITIALIZE
  CALL READ_FILE

  !CALL SMOOTH_COARSE_GRAIN
  CALL MAKE_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM OHcorr

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : ndata, coarse, nsteps, &
                         timestep, lifetime

  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: data_file, index_file

  CALL getarg(1, data_file)
  CALL getarg(2, index_file)

  OPEN(unit = 1, file = index_file)
  READ(1,*) ndata, nsteps, timestep
  CLOSE(1)

  ALLOCATE(coarse(nsteps,ndata)); coarse=0
  ALLOCATE(lifetime(nsteps)); lifetime = 0

  OPEN(unit = 1,file = data_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FILE
  USE histogram, ONLY : coarse, nsteps, ndata
  IMPLICIT NONE
  INTEGER :: step, i, frame
  
  DO step=1,nsteps

    READ(1,*) frame, ( coarse(step,i), i=1,ndata )

  END DO

END SUBROUTINE READ_FILE

SUBROUTINE MAKE_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : ndata, nsteps, coarse, lifetime 
  IMPLICIT NONE
  INTEGER                    :: frame1, frame2, idata
  INTEGER                    :: lt, c
  INTEGER, PARAMETER         :: nsep=100
  LOGICAL                    :: is_contiguous

  !frame1 is the duration of the interval 
  !DO frame1 = 0,nsteps/2,10
  DO frame1 = 0,nsteps/2

    lt=0
    c=0
!$omp parallel do reduction(+:lt,c)
    DO idata = 1,ndata

      !Change initial condition, look for combinations with same interval
      DO frame2 = 1,nsteps-frame1,nsep

        !IF ( is_contiguous( idata, frame2, frame2+frame1, 2 ) ) THEN 
        IF ( ALL( coarse(frame2:frame2+frame1,idata) == 1 ) ) THEN
        !IF ( ( coarse(frame2,idata) == 1 ) .and. ( coarse(frame1+frame2,idata) == 1 ) ) THEN
       
          lt = lt + 1
        
        ENDIF

        c = c + 1

      ENDDO

    ENDDO
!$omp end parallel do

    lifetime(frame1+1)   = dble(lt) / dble(c)

  ENDDO

END SUBROUTINE MAKE_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : lifetime, nsteps, timestep, ndata
  IMPLICIT NONE
  INTEGER                    :: i, j

  OPEN( unit = 2, file = "lifetime.dat" )

  DO i = 0,nsteps

    write(2,fmt = '(F15.7,3X,F15.7)'), &
         float(i)*timestep, lifetime(i+1) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE SMOOTH_COARSE_GRAIN
  !Avoid fast jumps in the coarse grain function
  USE histogram, ONLY : coarse, ndata, nsteps
  IMPLICIT NONE
  INTEGER                    :: idata, frame, initial, tout
  INTEGER, PARAMETER         :: forgiv=2 !10 
  
  DO idata = 1, ndata
  
    frame = 1
    tout  = 0
    initial = coarse(frame,idata)

    DO WHILE ( frame <= nsteps ) 

      IF ( coarse(frame,idata) /= initial ) THEN
        tout = tout + 1
      ELSEIF ( tout > 0 ) THEN 
        coarse(frame-tout:frame,idata) = initial
        tout = 0
      END IF
  
      !If the particle really changed layer
      IF ( tout > forgiv ) THEN

        initial = coarse(frame,idata)
        tout = 0

      END IF

      frame = frame + 1

    END DO

  END DO

END SUBROUTINE SMOOTH_COARSE_GRAIN

LOGICAL FUNCTION is_contiguous( idata, ti, tf, sep )
  !Test if atom is inside a layer from time ti to tf
  USE histogram, ONLY : coarse
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: idata, ti, tf, sep
  INTEGER                    :: initial

  initial = coarse(ti,idata)

  is_contiguous = ALL( coarse(ti:tf:sep,idata) == initial )

END FUNCTION is_contiguous
