MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
  DOUBLE PRECISION                         :: timestep
  INTEGER, ALLOCATABLE                     :: coarse(:,:)
  INTEGER                                  :: ndata, nsteps
END MODULE histogram

PROGRAM OHcorr
  IMPLICIT NONE

  CALL INITIALIZE
  CALL READ_FILE

  CALL SMOOTH_COARSE_GRAIN
  CALL MAKE_AUTOCORR 
  CALL AVERAGE_AUTOCORR
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM OHcorr

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : autocorr, ndata, coarse, nsteps, &
                         timestep

  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: data_file, index_file

  CALL getarg(1, data_file)
  CALL getarg(2, index_file)

  OPEN(unit = 1, file = index_file)
  READ(1,*) ndata, nsteps, timestep
  CLOSE(1)

  ALLOCATE(coarse(ndata,nsteps)); coarse=0
  ALLOCATE(autocorr(ndata,nsteps)); autocorr = 0

  OPEN(unit = 1,file = data_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FILE
  USE histogram, ONLY : coarse, nsteps, ndata
  IMPLICIT NONE
  INTEGER :: step, i, frame
  
  DO step=1,nsteps

    READ(1,*) frame, ( coarse(i,step), i=1,ndata )

  END DO

END SUBROUTINE READ_FILE

SUBROUTINE MAKE_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : ndata, nsteps, coarse, autocorr
  IMPLICIT NONE
  INTEGER                    :: frame1, frame2, idata
  INTEGER                    :: prod(ndata) 
  INTEGER, PARAMETER         :: nsep=20
  LOGICAL                    :: is_contiguous

  autocorr=0

  !frame1 is the duration of the interval 
  DO frame1 = 0,nsteps/2

    prod=0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nsteps-frame1,nsep

      DO idata = 1,ndata

        prod(idata) = prod(idata) + &
                    & coarse(idata, frame2) * &
                    & coarse(idata, frame2+frame1)

      ENDDO

    ENDDO

    autocorr(:,frame1+1) = dble(prod) / dble( frame2-1 ) * dble(nsep)

  ENDDO

END SUBROUTINE MAKE_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : autocorr, nsteps, timestep, ndata
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION           :: mean_autocorr(nsteps)

  OPEN( unit = 2, file = "autocorr.dat" )

  DO i = 0,nsteps-1

    !write(2,fmt = '(E11.4,3X,*(E14.7,3X))'), dble(i)*timestep, &
    !     ( autocorr(j,i+1) , j=1,ndata) 
    !     ( autocorr(j,i+1) / autocorr(j,1) , j=1,ndata) 
    write(2,fmt = '(E11.4,3X,E14.7)'), &
         dble(i)*timestep, autocorr(1,i+1) !/ autocorr(1,1)  

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE SMOOTH_COARSE_GRAIN
  !Avoid fast jumps in the coarse grain function
  USE histogram, ONLY : coarse, ndata, nsteps
  IMPLICIT NONE
  INTEGER                    :: idata, frame, initial, tout
  INTEGER, PARAMETER         :: forgiv=2
  
  DO idata = 1, ndata
  
    frame = 1
    tout  = 0
    initial = coarse(idata,frame)

    DO WHILE ( frame <= nsteps ) 

      IF ( coarse(idata,frame) /= initial ) THEN
        tout = tout + 1
      ELSEIF ( tout > 0 ) THEN 
        coarse(idata,frame-tout:frame) = initial
        tout = 0
      END IF
  
      !If the particle really changed layer
      IF ( tout > forgiv ) THEN

        initial = coarse(idata,frame)
        tout = 0

      END IF

      frame = frame + 1

    END DO

  END DO

END SUBROUTINE SMOOTH_COARSE_GRAIN

SUBROUTINE AVERAGE_AUTOCORR
  !Make first element of autocorr the average of
  ! normalized autocorrelation functions
  USE histogram, ONLY : autocorr, ndata, nsteps
  IMPLICIT NONE
  INTEGER :: idata, step, counter
  DOUBLE PRECISION :: buffer(nsteps)

  buffer=0.d0
  counter=0

  DO idata=1,ndata
    IF (autocorr(idata,1).ne.0.d0) then
      DO step = 1,nsteps
        autocorr(idata,step) = autocorr(idata,step) !/ autocorr(idata,1)
        buffer(step)=buffer(step)+autocorr(idata,step)
        counter = counter + 1
      END DO
    END IF
  END DO
      
  DO step = 1,nsteps
    autocorr(1,step) = buffer(step)/dble(counter)
  END DO

END SUBROUTINE AVERAGE_AUTOCORR

LOGICAL FUNCTION is_contiguous( idata, ti, tf, sep )
  !Test if atom is inside a layer from time ti to tf
  USE histogram, ONLY : coarse
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: idata, ti, tf
  INTEGER                    :: ts, initial
  INTEGER, OPTIONAL          :: sep

  ts = 1
  IF (present(sep)) ts=sep

  initial = coarse(idata,ti)

  is_contiguous = ALL( coarse(idata,ti:tf:ts) == initial )

END FUNCTION is_contiguous
