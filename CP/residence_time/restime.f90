!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
  INTEGER,PARAMETER                        :: nbins = 20
  DOUBLE PRECISION,PARAMETER               :: z0 = 9.00 , zf = 27.90 

END MODULE histogram

PROGRAM ResTime
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_POS 
    !CALL COARSE_GRAIN_EQUAL (frame)
    CALL COARSE_GRAIN_POS (frame)
  END DO
 
  CALL SMOOTH_COARSE_GRAIN
  CALL AUTOCORR_COARSE 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ResTime

SUBROUTINE INITIALIZE
  USE histogram , ONLY : autocorr, nbins 
  USE parameters, ONLY : natoms, nframes, nwater, coarse, pos
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  CALL READ_INDEX ( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(coarse(nframes,nwater)); coarse = 0
  ALLOCATE(autocorr(nframes,nbins)); autocorr = 0

  OPEN(unit = 1,file = file_name)
  
END SUBROUTINE INITIALIZE

SUBROUTINE AUTOCORR_COARSE
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : nbins, autocorr 
  USE parameters, ONLY : coarse, nframes, nwater
  IMPLICIT NONE
  INTEGER                    :: sum_coarse(nbins), sum_coarse_temp(nbins)
  INTEGER                    :: frame1, frame2, i_at, i_coord
  INTEGER                    :: layer, i_bin
  INTEGER                    :: counter(nbins), counter_tmp(nbins) 
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,nframes-1

    counter = 0
    sum_coarse = 0.D0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1

      counter_tmp = 0
      sum_coarse_temp = 0.D0

!$omp parallel do private(layer) reduction(+:counter_tmp, sum_dip_temp)
      DO i_at = 1,nwater

        layer = coarse(frame2,i_at)

        counter_tmp(layer) = counter_tmp(layer) + 1
        IF ( layer == coarse(frame2+frame1, i_at) ) &
           sum_coarse_temp(layer) = sum_coarse_temp(layer) + 1 

      ENDDO
!$omp end parallel do

      counter = counter + counter_tmp
      sum_coarse = sum_coarse + sum_coarse_temp

    ENDDO

    DO i_bin = 1, nbins

      IF ( counter(i_bin) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(frame1+1,i_bin) = dble(sum_coarse(i_bin)) / dble(counter(i_bin))
      END IF

    END DO

  ENDDO

END SUBROUTINE AUTOCORR_COARSE 

SUBROUTINE COARSE_GRAIN_EQUAL (frame)
  USE histogram,  ONLY : zf, z0, nbins
  USE parameters, ONLY : coarse, natoms, nwater, pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iw, iat, ind
  LOGICAL                    :: is_water_oxygen
  DOUBLE PRECISION           :: z

  iw = 1
  DO iat = 1, natoms

    IF ( is_water_oxygen(iat) ) THEN

      z = pos(iat,3) - nint(pos(iat,3)/box(3,3))*box(3,3)
      z = z + box(3,3)/2.d0 - z0

      ind = int(dble(nbins)*z/(zf-z0)) + 1
      coarse(frame,iw) = ind

      iw = iw +1

    END IF

  END DO

END SUBROUTINE COARSE_GRAIN_EQUAL

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : autocorr, nbins
  USE parameters, ONLY : nframes, dt, iprint 
  IMPLICIT NONE
  INTEGER                    :: layer, frame 
  DOUBLE PRECISION           :: z
  DOUBLE PRECISION,PARAMETER :: au_to_fs = 2.418884326509D-2
  DOUBLE PRECISION           :: tc  

  tc = au_to_fs*dt

  OPEN(unit = 2,file = 'restime.dat')

  DO frame = 1, nframes
  
    write(2,fmt = '(E11.4,*(3X,E12.5))'), frame*tc, autocorr(frame,:) / autocorr(1,:) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
