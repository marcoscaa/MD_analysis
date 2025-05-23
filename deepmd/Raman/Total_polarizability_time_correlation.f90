
MODULE correlation 
  IMPLICIT NONE 
  INTEGER                                  :: tcorr
  DOUBLE PRECISION, ALLOCATABLE            :: isopolar(:), anisopolar(:,:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: iso_corr(:)
  DOUBLE PRECISION, ALLOCATABLE            :: aniso_corr(:) 
END MODULE correlation

PROGRAM PTC 
  USE parameters, ONLY : nframes, coarse
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  DO frame = 1,nframes
    CALL READ_RAW_POLARIZABILITY 
    CALL COMPUTE_POLARIZABILITIES (frame) 
  END DO

  CALL POLAR_CORR 
  CALL SUBTRACT_CONSTANT
  CALL PRINT_RESULTS

  CLOSE(5)

END PROGRAM PTC

SUBROUTINE INITIALIZE
  USE correlation, ONLY  : isopolar, anisopolar, tcorr, &
                         iso_corr, aniso_corr 
  USE parameters, ONLY : polar, nwater, &
                         nframes, nequil, &
                         dt, mean_pol
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) nframes, nequil, dt, tcorr
  READ(1,*) mean_pol
  CLOSE(1)

  nwater=1

  ALLOCATE(polar(3,3,1)) ! polarizability tensor
  ALLOCATE(isopolar(nframes)) ! isotropic polarizability 
  ALLOCATE(anisopolar(3,3,nframes)) ! anisotropic polarizability 
  ALLOCATE(iso_corr(tcorr)); iso_corr=0.d0 
  ALLOCATE(aniso_corr(tcorr)); aniso_corr=0.d0

  OPEN(unit = 5,file = 'total_polarizability.raw')
  
END SUBROUTINE INITIALIZE

SUBROUTINE POLAR_CORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE correlation, ONLY : isopolar, anisopolar, iso_corr,  & 
                          aniso_corr,tcorr
  USE parameters,  ONLY : nframes, nwater
  IMPLICIT NONE
  DOUBLE PRECISION           :: sumiso
  DOUBLE PRECISION           :: sumaniso
  DOUBLE PRECISION           :: tracematmul
  INTEGER, PARAMETER         :: sep=1
  INTEGER                    :: frame1, frame2, i_at, j_at, nt0

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    nt0=0
    sumiso=0.;sumaniso=0.

    !Change initial condition, look for combinations with same interval
!$omp parallel do reduction(+:sumiso,sumaniso,nt0)
    DO frame2 = 1,nframes-frame1,sep

      nt0=nt0+1

      sumiso = sumiso + &
          isopolar(frame2) * isopolar(frame2+frame1) 

      sumaniso = sumaniso + &
          tracematmul(anisopolar(:,:,frame2), anisopolar(:,:,frame2+frame1))

    ENDDO
!$omp end parallel do

    !Time-averaged correlation functions
    iso_corr(frame1+1)   = sumiso   / dble(nt0)
    aniso_corr(frame1+1) = sumaniso / dble(nt0)

  ENDDO

END SUBROUTINE POLAR_CORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE correlation, ONLY : iso_corr, aniso_corr, tcorr
  USE parameters,  ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN(unit = 10,file = 'ptc.dat')

  WRITE(10, *), "# Time, Iso Intra-corr, Iso Inter-corr, Aniso Intra-corr, Aniso Inter-corr"

  DO i = 1,tcorr

    write(10,fmt = '(F12.8,2(3X,E18.11))'), dble(i-1)*dt, &
          iso_corr(i), aniso_corr(i) 

  ENDDO

  CLOSE(10)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE COMPUTE_POLARIZABILITIES (frame)
  !Compute isotropic and anisotropic polarizabilities
  USE correlation, ONLY : isopolar, anisopolar
  USE parameters,  ONLY : polar, nwater, mean_pol
  IMPLICIT NONE
  INTEGER, INTENT(in)        :: frame
  INTEGER                    :: iw, ii

  isopolar(frame) = ( polar(1,1,1) + polar(2,2,1) + polar(3,3,1) ) / 3.
  anisopolar(:,:,frame) = polar(:,:,1)

  DO ii=1,3
    anisopolar(ii,ii,frame) = anisopolar(ii,ii,frame) - isopolar(frame)
  END DO

END SUBROUTINE COMPUTE_POLARIZABILITIES

SUBROUTINE SUBTRACT_CONSTANT
  USE correlation, ONLY : isopolar, anisopolar, & 
                          iso_corr, aniso_corr, tcorr
  USE parameters,  ONLY : nframes, nwater
  IMPLICIT NONE 
  INTEGER                    :: it, iw
  DOUBLE PRECISION           :: mean_iso, mean_aniso(3,3)
  DOUBLE PRECISION           :: tracematmul
   
  mean_iso=0.d0;mean_aniso=0.d0

  DO it=1,nframes

      mean_iso = mean_iso + isopolar(it)
      mean_aniso = mean_aniso + anisopolar(:,:,it)

  END DO

  mean_iso = mean_iso / dble(nframes)
  mean_aniso = mean_aniso / dble(nframes)

  iso_corr   = iso_corr &
             - mean_iso * mean_iso
  aniso_corr = aniso_corr &
             - tracematmul(mean_aniso,mean_aniso)

END SUBROUTINE SUBTRACT_CONSTANT

DOUBLE PRECISION FUNCTION tracematmul( A, B )
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(in) :: A(3,3), B(3,3)

  tracematmul = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1) &
              + A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2) &
              + A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3) 

END FUNCTION tracematmul
