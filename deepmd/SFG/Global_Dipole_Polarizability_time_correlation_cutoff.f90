
MODULE correlation 
  IMPLICIT NONE 
  INTEGER                         :: tcorr !Length of the time corr. function
  DOUBLE PRECISION, ALLOCATABLE   :: polart(:,:,:), dipolet(:,:)
  DOUBLE PRECISION, ALLOCATABLE   :: corr(:,:,:)
END MODULE correlation

PROGRAM PTC 
  USE parameters, ONLY : nframes, coarse
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  !CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_TOTAL_DIPOLE_POLARIZABILITY_RAW(frame)
  END DO

  CALL DIPOLE_POLAR_CORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM PTC

SUBROUTINE INITIALIZE
  USE correlation, ONLY : polart, dipolet, tcorr, &
                          corr
  USE parameters,  ONLY : nframes, nequil, dt
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) nframes, nequil, dt, tcorr
  CLOSE(1)

  ALLOCATE(polart(2,2,nframes)) ! total polarizability 
  ALLOCATE(dipolet(2,nframes)) ! total dipole
  ALLOCATE(corr(2,2,tcorr)); corr=0.d0 

  CLOSE(1)

  OPEN(unit = 1,file = 'total_dipole.raw')
  OPEN(unit = 2,file = 'total_polarizability.raw')
  
END SUBROUTINE INITIALIZE

SUBROUTINE DIPOLE_POLAR_CORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE correlation, ONLY : polart, dipolet, tcorr, & 
                          corr
  USE parameters,  ONLY : nframes
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_up1, sum_down1
  DOUBLE PRECISION           :: sum_up2, sum_down2
  INTEGER, PARAMETER         :: sep=1
  INTEGER                    :: frame1, frame2, nt0
  INTEGER                    :: frame1p2

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    nt0=0
    sum_up1=0.d0;sum_down1=0.d0
    sum_up2=0.d0;sum_down2=0.d0

    !Change initial condition, look for combinations with same interval
!$omp parallel do private(frame1p2) reduction(+:sum_up1,sum_up2,sum_down1,sum_down2,nt0)
    DO frame2 = 1,nframes-frame1,sep

      frame1p2=frame1+frame2
      nt0=nt0+1

      sum_up1 = sum_up1 + &
          dipolet(1,frame2) * polart(1,1,frame1p2) 
      sum_up2 = sum_up2 + &
          dipolet(1,frame2) * polart(1,2,frame1p2) 
      sum_down1 = sum_down1 + &
          dipolet(2,frame2) * polart(2,1,frame1p2) 
      sum_down2 = sum_down2 + &
          dipolet(2,frame2) * polart(2,2,frame1p2) 

    ENDDO
!$omp end parallel do

    !Time-averaged correlation functions
    corr(1,1,frame1+1) = sum_up1 / dble(nt0)
    corr(1,2,frame1+1) = sum_up2 / dble(nt0)
    corr(2,1,frame1+1) = sum_down1 / dble(nt0)
    corr(2,2,frame1+1) = sum_down2 / dble(nt0)

  ENDDO

END SUBROUTINE DIPOLE_POLAR_CORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE correlation, ONLY : corr, tcorr
  USE parameters,  ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN(unit = 10,file = 'total_dptc.dat')

  WRITE(10, *), "# Time, Corr_xxz(layer1), Corr_yyz(layer1),Corr_xxz(layer2), Corr_yyz(layer2)"

  DO i = 1,tcorr

    write(10,fmt = '(F12.8,4(3X,E18.11))'), dble(i-1)*dt, &
          corr(1,1,i), corr(1,2,i), corr(2,1,i), corr(2,2,i)

  ENDDO

  CLOSE(10)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE READ_TOTAL_DIPOLE_POLARIZABILITY_RAW (frame)
  USE correlation, ONLY : dipolet, polart
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: i

  READ(1,*) dipolet(1,frame), dipolet(2,frame)
  READ(2,*) polart(1,1,frame), polart(1,2,frame), & 
            polart(2,1,frame), polart(2,2,frame)

END SUBROUTINE READ_TOTAL_DIPOLE_POLARIZABILITY_RAW
