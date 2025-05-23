
MODULE correlation 
  IMPLICIT NONE 
  INTEGER                                  :: tcorr
  DOUBLE PRECISION, ALLOCATABLE            :: isopolar(:,:), anisopolar(:,:,:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: iso_intra_corr(:), iso_inter_corr(:) 
  DOUBLE PRECISION, ALLOCATABLE            :: aniso_intra_corr(:), aniso_inter_corr(:) 
  DOUBLE PRECISION, PARAMETER              :: rcut=6.d0
END MODULE correlation

PROGRAM PTC 
  USE parameters, ONLY : nframes, coarse
  USE correlation, ONLY: rcut
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_RAW_POLARIZABILITY 
    CALL READ_RAW_POS_BOX
    CALL COMPUTE_POLARIZABILITIES (frame) 
    CALL MAKE_NEIGHBOR_LIST (frame,rcut)
  END DO

  CALL SUBTRACT_AVERAGE
  CALL POLAR_CORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM PTC

SUBROUTINE INITIALIZE
  USE correlation, ONLY  : isopolar, anisopolar, tcorr, &
                         iso_intra_corr, iso_inter_corr, &
                         aniso_intra_corr, aniso_inter_corr 
  USE parameters, ONLY : pos, polar, natoms, nwater, &
                         nframes, nequil, within_cutoff, &
                         dt, atype, mean_pol
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nwater, dt, tcorr
  READ(1,*) mean_pol
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(polar(3,3,nwater)) ! molecular polarizability tensors
  ALLOCATE(isopolar(nwater,nframes)) ! isotropic molecular polarizability 
  ALLOCATE(anisopolar(3,3,nwater,nframes)) ! anisotropic molecular polarizability 
  ALLOCATE(within_cutoff(nwater,nwater,nframes))
  ALLOCATE(iso_intra_corr(tcorr)); iso_intra_corr=0.d0 
  ALLOCATE(iso_inter_corr(tcorr)); iso_inter_corr=0.d0
  ALLOCATE(aniso_intra_corr(tcorr)); aniso_intra_corr=0.d0
  ALLOCATE(aniso_inter_corr(tcorr)); aniso_inter_corr=0.d0
  ALLOCATE(atype(natoms))

  !Read atype from type.raw
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  CLOSE(1)

  OPEN(unit = 1,file = 'coord.raw')
  OPEN(unit = 2,file = 'box.raw')
  OPEN(unit = 5,file = 'polarizability.raw')
  
END SUBROUTINE INITIALIZE

SUBROUTINE POLAR_CORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE correlation, ONLY : isopolar, anisopolar, iso_intra_corr, iso_inter_corr, & 
                         aniso_intra_corr, aniso_inter_corr, tcorr
  USE parameters,  ONLY : nframes, nwater, within_cutoff
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

          IF ( within_cutoff(j_at, i_at, frame2) ) THEN

            sumiso_inter = sumiso_inter + &
                isopolar(i_at,frame2) * isopolar(j_at,frame2+frame1) + &
                isopolar(j_at,frame2) * isopolar(i_at,frame2+frame1) 

            sumaniso_inter = sumaniso_inter + &
                tracematmul(anisopolar(:,:,i_at,frame2), anisopolar(:,:,j_at,frame2+frame1)) + &
                tracematmul(anisopolar(:,:,j_at,frame2), anisopolar(:,:,i_at,frame2+frame1)) 

          END IF

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

  OPEN(unit = 5,file = 'ptc.dat')

  WRITE(5, *), "# Time, Iso Intra-corr, Iso Inter-corr, Aniso Intra-corr, Aniso Inter-corr"

  DO i = 1,tcorr

    write(5,fmt = '(F12.8,4(3X,E18.11))'), dble(i-1)*dt, &
          iso_intra_corr(i), iso_inter_corr(i), &
          aniso_intra_corr(i), aniso_inter_corr(i) 

  ENDDO

  CLOSE(5)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE SUBTRACT_AVERAGE
  USE correlation, ONLY : isopolar, anisopolar
  USE parameters,  ONLY : nframes, nwater
  IMPLICIT NONE 
  INTEGER                    :: it, iw
  DOUBLE PRECISION           :: mean_iso, mean_aniso(3,3)
   
  mean_iso=0.d0;mean_aniso=0.d0

  DO it=1,nframes
    DO iw=1,nwater

      mean_iso = mean_iso + isopolar(iw,it)
      mean_aniso = mean_aniso + anisopolar(:,:,iw,it)

    END DO
  END DO

  isopolar = isopolar - mean_iso / dble(nwater*nframes)
  mean_aniso = mean_aniso / dble(nwater*nframes)

  DO it=1,nframes
    DO iw=1,nwater
      anisopolar(:,:,iw,it) = anisopolar(:,:,iw,it) - mean_aniso 
    END DO
  END DO

END SUBROUTINE SUBTRACT_AVERAGE

SUBROUTINE COMPUTE_POLARIZABILITIES (frame)
  !Compute isotropic and anisotropic polarizabilities
  USE correlation, ONLY : isopolar, anisopolar
  USE parameters,  ONLY : polar, nwater
  IMPLICIT NONE
  INTEGER, INTENT(in)        :: frame
  INTEGER                    :: iw, i, j

  DO iw=1,nwater
    
    isopolar(iw,frame) = ( polar(1,1,iw) + polar(2,2,iw) + polar(3,3,iw) ) / 3.

    DO i=1,3
      DO j=i+1,3
        anisopolar(i,j,iw,frame) = ( polar(i,j,iw) + polar(j,i,iw) ) / 2.
        anisopolar(j,i,iw,frame) = anisopolar(i,j,iw,frame) 
      ENDDO
        anisopolar(i,i,iw,frame) = polar(i,i,iw) - isopolar(iw,frame)
    ENDDO

  END DO

END SUBROUTINE COMPUTE_POLARIZABILITIES

DOUBLE PRECISION FUNCTION tracematmul( A, B )
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(in) :: A(3,3), B(3,3)

  tracematmul = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1) &
              + A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2) &
              + A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3) 

END FUNCTION tracematmul
