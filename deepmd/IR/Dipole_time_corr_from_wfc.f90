
MODULE histogram 
  IMPLICIT NONE 
  INTEGER                                  :: tcorr
  DOUBLE PRECISION, ALLOCATABLE            :: dipt(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr_inter(:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr_intra(:)
  DOUBLE PRECISION, PARAMETER              :: rcut=6.d0
END MODULE histogram

PROGRAM DAC 
  USE parameters, ONLY : nframes, coarse
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_RAW_WANNIER 
    CALL READ_RAW_POS_BOX
    CALL ASSIGN_DIPOLE (frame) 
    CALL MAKE_NEIGHBOR_LIST (frame)
  END DO

  CALL DIPOLE_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM DAC

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : dipt, autocorr_inter, autocorr_intra, tcorr 
  USE parameters, ONLY : pos, wannier, natoms, nwater, &
                         nframes, nequil, within_cutoff, dt, &
                         atype
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nwater, dt, tcorr
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(wannier(3,4,nwater)) 
  ALLOCATE(dipt(3,nwater,nframes)) ! dipole moment of all water
  ALLOCATE(within_cutoff(nwater,nwater,nframes))
  ALLOCATE(autocorr_inter(nframes)); autocorr_inter = 0
  ALLOCATE(autocorr_intra(nframes)); autocorr_intra = 0
  ALLOCATE(atype(natoms))

  !Read atype from type.raw
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  CLOSE(1)

  OPEN(unit = 1,file = 'coord.raw')
  OPEN(unit = 2,file = 'box.raw')
  OPEN(unit = 4,file = 'wannier.raw')
  
END SUBROUTINE INITIALIZE

SUBROUTINE ASSIGN_DIPOLE(frame)
  USE histogram, ONLY : dipt
  USE parameters, ONLY : pos, wannier, natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, iw, iwfc, iO
  DOUBLE PRECISION                     :: d(4)
  LOGICAL                              :: is_water_oxygen
  
  iw = 0
  dipt(:,:,frame) = 0

  DO iat=1,natoms

    !Electronic part
    IF ( is_water_oxygen(iat) ) THEN
  
      iw = iw + 1
      DO iwfc = 1, 4
        !Assuming Wannier coordinates are relative to oxygen atom
        dipt(:,iw,frame) = dipt(:,iw,frame) - 2.d0*wannier(:,iwfc,iw)
      END DO
  
      iO = iat

    ELSE

      CALL DISTANCE_VECTOR( pos(:,iO), pos(:,iat), d ) 
      !Ionic part
      dipt(:,iw,frame) = dipt(:,iw,frame) + d(1:3) 

    END IF

  END DO

END SUBROUTINE ASSIGN_DIPOLE

SUBROUTINE DIPOLE_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : dipt, autocorr_inter, autocorr_intra, tcorr, rcut
  USE parameters, ONLY : nframes, nwater, within_cutoff
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_intra, sum_inter
  DOUBLE PRECISION           :: dotprod
  INTEGER, PARAMETER         :: sep=10
  INTEGER                    :: frame1, frame2, i_at, j_at, n0

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    n0=0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,sep

      sum_inter = 0.D0
      sum_intra = 0.D0
      n0=n0+1

!$omp parallel do private(dotprod) reduction(+:sum_inter,sum_intra)
      DO i_at = 1,nwater

        sum_intra = sum_intra + &
                    sum( dipt(:,i_at,frame2) * dipt(:,i_at, frame2+frame1) )

        DO j_at = i_at+1, nwater

          !IF ( within_cutoff(i_at, j_at, frame2) ) THEN
          !IF ( Dist_time(i_at, j_at, frame2,frame2+frame1) < rcut ) THEN

            dotprod = sum( dipt(:,i_at,frame2) * &
                           dipt(:,j_at, frame2+frame1) ) + & 
                      sum( dipt(:,j_at,frame2) * &
                           dipt(:,i_at, frame2+frame1) ) 

            sum_inter = sum_inter + dotprod

          !END IF

        ENDDO

      ENDDO
!$omp end parallel do
 
      autocorr_inter(frame1+1) = autocorr_inter(frame1+1) + sum_inter
      autocorr_intra(frame1+1) = autocorr_intra(frame1+1) + sum_intra

    ENDDO

    !Time-averaged autocorrelation function
    autocorr_inter(frame1+1) = autocorr_inter(frame1+1) / dble(n0)
    autocorr_intra(frame1+1) = autocorr_intra(frame1+1) / dble(n0)

  ENDDO

END SUBROUTINE DIPOLE_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,   ONLY : autocorr_inter, autocorr_intra
  USE parameters, ONLY : box, dt, nlayers, nframes
  IMPLICIT NONE
  INTEGER                    :: i, j

  OPEN(unit = 5,file = 'dac.dat')

  DO i = 0,nframes-1

    write(5,fmt = '(E18.11,2(3X,E12.5))'), dble(i)*dt, autocorr_inter(i+1), &
                                                       autocorr_intra(i+1) 

  ENDDO

  CLOSE(5)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE MAKE_NEIGHBOR_LIST(frame)
  USE PARAMETERS, ONLY : natoms, nwater, within_cutoff
  USE histogram, ONLY : rcut
  IMPLICIT NONE
  INTEGER, INTENT(in)        :: frame
  INTEGER                    :: iat, jat, iw, jw
  REAL*8                     :: d, Dist
  LOGICAL                    :: is_water_oxygen
  
  iw = 0
  jw=0

  DO iat = 1,natoms 
 
    IF (is_water_oxygen(iat)) THEN
      
      iw = iw + 1

      DO jat = iat,natoms

        IF (is_water_oxygen(jat)) THEN

          jw = jw + 1

          IF (iat==jat) THEN 
            within_cutoff(iw,jw,frame) = .true.
          ELSE
            IF ( Dist(iat,jat) < rcut ) THEN
              within_cutoff(iw,jw,frame) = .true.
            ELSE
              within_cutoff(iw,jw,frame) = .false.
            ENDIF
            within_cutoff(iw,jw,frame) = within_cutoff(jw,iw,frame)
          ENDIF
              
        END IF

      END DO

    END IF
  
  END DO

END SUBROUTINE MAKE_NEIGHBOR_LIST
