
MODULE histogram 
  IMPLICIT NONE 
  INTEGER                                  :: tcorr
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr_inter(:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr_intra(:)
  DOUBLE PRECISION, ALLOCATABLE            :: dipolet(:,:,:)
  DOUBLE PRECISION, PARAMETER              :: rcut=6.d0
END MODULE histogram

PROGRAM DAC 
  USE parameters, ONLY : nframes, coarse
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_RAW_DIPOLE (frame) 
    CALL READ_RAW_POS_BOX
    CALL ASSIGN_DIPOLE (frame) 
    !CALL MAKE_NEIGHBOR_LIST (frame)
  END DO

  CALL DIPOLE_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM DAC

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : autocorr_inter, autocorr_intra, tcorr, &
                         dipolet 
  USE parameters, ONLY : pos, natoms, nwater, nframes, nequil, &
                         within_cutoff, dt, atype, dipole
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nwater, dt, tcorr
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(dipole(3,nwater)) ! dipole moment of all water
  ALLOCATE(dipolet(3,nwater,nframes)) ! dipole moment of all water
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
  OPEN(unit = 4,file = 'dipole.raw')
  
END SUBROUTINE INITIALIZE

SUBROUTINE ASSIGN_DIPOLE(frame)
  !We assume atomic type to be on the order: O H H O H H
  USE parameters, ONLY : dipole, pos, nwater, atype, natoms
  USE histogram , ONLY : dipolet
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, ih, iw, ind_H(2)
  DOUBLE PRECISION                     :: d(4), diptmp(3)
  
  iw = 0

  DO iat=1,natoms

    IF (atype(iat)==0) THEN

      iw = iw + 1 

      !Electronic part
      diptmp = -2.d0*dipole(:,iw)
      !dipole(:,iw,frame) = 0.d0 

      CALL get_H2(iat, ind_H)

      DO ih=1,2 
        
        CALL DISTANCE_VECTOR(pos(:,iat), pos(:,ind_H(ih)), d)
        diptmp = diptmp + d(1:3) 

      END DO

      CALL APPLY_PBC(diptmp, dipolet(:,iw,frame))

    END IF

  END DO

END SUBROUTINE ASSIGN_DIPOLE

SUBROUTINE DIPOLE_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : autocorr_inter, autocorr_intra, tcorr, rcut, dipolet
  USE parameters, ONLY : nframes, nwater, within_cutoff
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_intra, sum_inter
  DOUBLE PRECISION           :: dotprod
  INTEGER, PARAMETER         :: sep=10
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

  DO i = 0,tcorr-1

    write(5,fmt = '(E18.11,2(3X,E18.11))'), dble(i)*dt, autocorr_inter(i+1), &
                                                       autocorr_intra(i+1) 

  ENDDO

  CLOSE(5)

END SUBROUTINE PRINT_RESULTS
