
MODULE histogram 
  IMPLICIT NONE 
  INTEGER                                  :: tcorr
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr_inter(:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr_intra(:)
  DOUBLE PRECISION, ALLOCATABLE            :: dipolet(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: dipole(:,:)
  DOUBLE PRECISION, PARAMETER              :: rcut=6.d0
END MODULE histogram

PROGRAM DAC 
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_DIPOLE 
    CALL ASSIGN_DIPOLE(frame)
  END DO

  CALL DIPOLE_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM DAC

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : autocorr_inter, autocorr_intra, tcorr, &
                         dipolet, dipole 
  USE parameters, ONLY : natoms, nframes, nequil, &
                         dt, atype, nwater
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file, dipole_file

  CALL getarg(1, dipole_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nwater, nframes, nequil, dt, tcorr
  CLOSE(1)

  ALLOCATE(dipole(3,natoms)) ! dipole moment of all water
  ALLOCATE(dipolet(3,nwater,nframes)) ! dipole moment of all water
  !ALLOCATE(within_cutoff(natoms,natoms,nframes))
  ALLOCATE(autocorr_inter(nframes)); autocorr_inter = 0
  ALLOCATE(autocorr_intra(nframes)); autocorr_intra = 0

  OPEN(unit = 1,file = dipole_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_DIPOLE 
  !Read LAMMPS atom file
  USE parameters, ONLY : natoms, box
  USE histogram, ONLY: dipole
  IMPLICIT NONE
  INTEGER                    :: iat, ind, atp
  REAL*8                     :: box_tmp(2), box_min(3)
 
  DO iat=1,5
    READ(1,*)
  END DO

  !Read Box
  DO iat=1,3
    READ(1,*) box_tmp(1), box_tmp(2)
    box(iat) = box_tmp(2)-box_tmp(1)
    box_min(iat)=box_tmp(1)
  END DO

  READ(1,*)

  DO iat = 1, natoms
    READ(1,*) ind, dipole(1,ind), dipole(2,ind), dipole(3,ind)
  END DO

END SUBROUTINE READ_DIPOLE

SUBROUTINE ASSIGN_DIPOLE(frame)
  USE parameters, ONLY : dipole=>pos, nwater
  USE histogram , ONLY : dipolet
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat
  
  DO iat = 1,nwater
    dipolet(:,iat,frame) = dipole(:,iat)
  END DO

END SUBROUTINE ASSIGN_DIPOLE

SUBROUTINE DIPOLE_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : autocorr_inter, autocorr_intra, tcorr, rcut, dipolet
  USE parameters, ONLY : nframes, nwater
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_intra, sum_inter
  DOUBLE PRECISION           :: dotprod
  INTEGER, PARAMETER         :: sep=100
  INTEGER                    :: frame1, frame2, i_at, j_at, n0, ipol

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    n0=0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,sep

      sum_inter = 0.D0
      sum_intra = 0.D0
      n0=n0+1

!$omp parallel do reduction(+:sum_inter,sum_intra)
      DO i_at = 1,nwater

        DO ipol=1,3
          sum_intra = sum_intra + &
                      dipolet(ipol,i_at,frame2) * dipolet(ipol,i_at, frame2+frame1)
        END DO

        DO j_at = i_at+1, nwater

          !IF ( within_cutoff(i_at, j_at, frame2) ) THEN
          !IF ( Dist_time(i_at, j_at, frame2,frame2+frame1) < rcut ) THEN

          DO ipol=1,3
            sum_inter = sum_inter + dipolet(ipol,i_at,frame2) * &
                                    dipolet(ipol,j_at, frame2+frame1) + & 
                                    dipolet(ipol,j_at,frame2) * &
                                    dipolet(ipol,i_at, frame2+frame1)  
          END DO

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
  USE histogram,   ONLY : autocorr_inter, autocorr_intra, tcorr
  USE parameters, ONLY : box, dt, nlayers, nframes
  IMPLICIT NONE
  INTEGER                    :: i, j

  OPEN(unit = 5,file = 'dac.dat')

  write(5, fmt='(27A)') '#Time DAC(inter) DAC(intra)'

  DO i = 0,tcorr-1

    write(5,fmt = '(E18.11,2(3X,E18.11))') dble(i)*dt, autocorr_inter(i+1), &
                                                       autocorr_intra(i+1) 

  ENDDO

  CLOSE(5)

END SUBROUTINE PRINT_RESULTS
