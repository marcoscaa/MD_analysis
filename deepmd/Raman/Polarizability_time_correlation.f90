
MODULE correlation 
  IMPLICIT NONE 
  INTEGER                                  :: tcorr
  DOUBLE PRECISION, ALLOCATABLE            :: isopolar(:,:), anisopolar(:,:,:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: iso_intra_corr(:), iso_inter_corr(:) 
  DOUBLE PRECISION, ALLOCATABLE            :: aniso_intra_corr(:), aniso_inter_corr(:) 
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
                         iso_intra_corr, iso_inter_corr, &
                         aniso_intra_corr, aniso_inter_corr 
  USE parameters, ONLY : polar, nwater, &
                         nframes, nequil, &
                         dt, mean_pol
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) nwater, nframes, nequil, dt, tcorr
  READ(1,*) mean_pol
  CLOSE(1)

  ALLOCATE(polar(3,3,nwater)) ! molecular polarizability tensors
  ALLOCATE(isopolar(nwater,nframes)) ! isotropic molecular polarizability 
  ALLOCATE(anisopolar(3,3,nwater,nframes)) ! anisotropic molecular polarizability 
  ALLOCATE(iso_intra_corr(tcorr)); iso_intra_corr=0.d0 
  ALLOCATE(iso_inter_corr(tcorr)); iso_inter_corr=0.d0
  ALLOCATE(aniso_intra_corr(tcorr)); aniso_intra_corr=0.d0
  ALLOCATE(aniso_inter_corr(tcorr)); aniso_inter_corr=0.d0

  OPEN(unit = 5,file = 'polarizability.raw')
  
END SUBROUTINE INITIALIZE

SUBROUTINE POLAR_CORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE correlation, ONLY : isopolar, anisopolar, iso_intra_corr, iso_inter_corr, & 
                         aniso_intra_corr, aniso_inter_corr, tcorr
  USE parameters,  ONLY : nframes, nwater
  IMPLICIT NONE
  DOUBLE PRECISION           :: sumiso_intra, sumiso_inter
  DOUBLE PRECISION           :: sumaniso_intra, sumaniso_inter
  DOUBLE PRECISION           :: tracematmul
  INTEGER, PARAMETER         :: sep=10
  INTEGER                    :: frame1, frame2, i_at, j_at, nt0

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    nt0=0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,sep

      sumiso_intra=0.d0;sumiso_inter=0.d0
      sumaniso_intra=0.d0;sumaniso_inter=0.d0
      nt0=nt0+1

!$omp parallel do reduction(+:sumiso_inter,sumiso_intra,sumaniso_intra,sumaniso_inter)
      DO i_at = 1,nwater

        sumiso_intra = sumiso_intra + &
            isopolar(i_at,frame2) * isopolar(i_at,frame2+frame1) 

        sumaniso_intra = sumaniso_intra + &
            tracematmul(anisopolar(:,:,i_at,frame2), anisopolar(:,:,i_at,frame2+frame1))

        DO j_at = i_at+1, nwater

          sumiso_inter = sumiso_inter + &
              isopolar(i_at,frame2) * isopolar(j_at,frame2+frame1) + &
              isopolar(j_at,frame2) * isopolar(i_at,frame2+frame1) 

          sumaniso_inter = sumaniso_inter + &
              tracematmul(anisopolar(:,:,i_at,frame2), anisopolar(:,:,j_at,frame2+frame1)) + &
              tracematmul(anisopolar(:,:,j_at,frame2), anisopolar(:,:,i_at,frame2+frame1)) 

        ENDDO

      ENDDO
!$omp end parallel do

      iso_intra_corr(frame1+1)   = iso_intra_corr(frame1+1)   + sumiso_intra
      iso_inter_corr(frame1+1)   = iso_inter_corr(frame1+1)   + sumiso_inter
      aniso_intra_corr(frame1+1) = aniso_intra_corr(frame1+1) + sumaniso_intra
      aniso_inter_corr(frame1+1) = aniso_inter_corr(frame1+1) + sumaniso_inter

    ENDDO

    !Time-averaged correlation functions
    iso_intra_corr(frame1+1)   = iso_intra_corr(frame1+1)   / dble(nt0)
    iso_inter_corr(frame1+1)   = iso_inter_corr(frame1+1)   / dble(nt0)
    aniso_intra_corr(frame1+1) = aniso_intra_corr(frame1+1) / dble(nt0)
    aniso_inter_corr(frame1+1) = aniso_inter_corr(frame1+1) / dble(nt0)

  ENDDO

END SUBROUTINE POLAR_CORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE correlation, ONLY : iso_intra_corr, iso_inter_corr, & 
                          aniso_intra_corr, aniso_inter_corr, &
                          tcorr
  USE parameters,  ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN(unit = 10,file = 'ptc.dat')

  WRITE(10, *), "# Time, Iso Intra-corr, Iso Inter-corr, Aniso Intra-corr, Aniso Inter-corr"

  DO i = 1,tcorr

    write(10,fmt = '(F12.8,4(3X,E18.11))'), dble(i-1)*dt, &
          iso_intra_corr(i), iso_inter_corr(i), &
          aniso_intra_corr(i), aniso_inter_corr(i) 

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

  DO iw=1,nwater
    
    isopolar(iw,frame) = ( polar(1,1,iw) + polar(2,2,iw) + polar(3,3,iw) ) / 3.
    anisopolar(:,:,iw,frame) = polar(:,:,iw)

    DO ii=1,3
      anisopolar(ii,ii,iw,frame) = anisopolar(ii,ii,iw,frame) - isopolar(iw,frame)
    END DO

    isopolar(iw,frame) = isopolar(iw,frame) !+ mean_pol

  END DO

END SUBROUTINE COMPUTE_POLARIZABILITIES

SUBROUTINE SUBTRACT_CONSTANT
  USE correlation, ONLY : isopolar, anisopolar, iso_intra_corr, iso_inter_corr, & 
                         aniso_intra_corr, aniso_inter_corr, tcorr
  USE parameters,  ONLY : nframes, nwater
  IMPLICIT NONE 
  INTEGER                    :: it, iw
  DOUBLE PRECISION           :: mean_iso, mean_aniso(3,3)
  DOUBLE PRECISION           :: tracematmul
   
  mean_iso=0.d0;mean_aniso=0.d0

  DO it=1,nframes
    DO iw=1,nwater

      mean_iso = mean_iso + isopolar(iw,it)
      mean_aniso = mean_aniso + anisopolar(:,:,iw,it)

    END DO
  END DO

  mean_iso = mean_iso / dble(nwater*nframes)
  mean_aniso = mean_aniso / dble(nwater*nframes)

  iso_intra_corr   = iso_intra_corr &
                   - dble(nwater) * mean_iso * mean_iso
  iso_inter_corr   = iso_inter_corr &
                   - dble(nwater*(nwater-1)) * mean_iso * mean_iso
  aniso_intra_corr = aniso_intra_corr &
                   - dble(nwater) * tracematmul(mean_aniso,mean_aniso)
  aniso_inter_corr = aniso_inter_corr &
                   - dble(nwater*(nwater-1)) * tracematmul(mean_aniso,mean_aniso)

END SUBROUTINE SUBTRACT_CONSTANT

DOUBLE PRECISION FUNCTION tracematmul( A, B )
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(in) :: A(3,3), B(3,3)

  tracematmul = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1) &
              + A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2) &
              + A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3) 

END FUNCTION tracematmul
