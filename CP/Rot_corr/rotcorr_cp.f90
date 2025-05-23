!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: dipole_dist(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)

END MODULE histogram 

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_POS
    CALL COARSE_GRAIN_POS (frame)
    CALL DIPOLE_DISTRIBUTION (frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL ROT_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : autocorr, dipole_dist
  USE parameters, ONLY : natoms, nframes, pos, &
                         nwater, nlayers, coarse
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  CALL READ_INDEX( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(coarse(nframes,nwater)); coarse=1
  ALLOCATE(autocorr(nframes,nlayers)); autocorr = 0.d0
  ALLOCATE(dipole_dist(nframes,nwater,3)); dipole_dist=0.d0

  OPEN(unit = 1,file = pos_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE DIPOLE_DISTRIBUTION (frame)
  !Generate the dipole distribution for each water in a specific frame
  USE histogram,  ONLY : dipole_dist
  USE parameters, ONLY : pos, natoms, ind_atom
  IMPLICIT NONE
  INTEGER                    :: frame, i, j
  INTEGER                    :: ind_O, ind_H(2)
  DOUBLE PRECISION           :: dipole
  LOGICAL                    :: is_water_oxygen

  !Water oxygen atom index
  ind_O = 0

  DO i = 1,natoms

    !Calculate the dipole as:
    !dipole = ( ( O - H1 ) + ( O - H2) ) / norm(dipole)
    IF ( is_water_oxygen(i) ) THEN

      ind_O = ind_O + 1
      ind_H = 0
      CALL get_H(i, ind_H)
   
      !Unit dipole vector taking PBC into account 
      DO j = 1,3 
        dipole_dist (frame,ind_O,j) = & 
         dipole(i,ind_H(1),ind_H(2),j)  
      END DO

    END IF

  ENDDO

END SUBROUTINE DIPOLE_DISTRIBUTION

SUBROUTINE ROT_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : dipole_dist, autocorr
  USE parameters, ONLY : nframes, nwater, coarse, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_dip(nlayers), surv_prob(nlayers) 
  DOUBLE PRECISION           :: second_leg, dotprod
  INTEGER                    :: frame1, frame2, i_at, i_coord
  INTEGER                    :: layer, i_bin
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers) 
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,nframes-1

    surv_prob = 0.d0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1

      counter_in = 0
      counter_all = 0
      sum_dip = 0.D0

!$omp parallel do private(layer,dotprod) reduction(+:counter_in, counter_all, sum_dip)
      DO i_at = 1,nwater

        layer = coarse(frame2,i_at)
        counter_all(layer) = counter_all(layer) + 1

        IF ( in_layer(i_at,frame2,frame2+frame1) ) THEN

          counter_in(layer) = counter_in(layer) + 1
          dotprod = 0.d0

          DO i_coord = 1,3
            dotprod = dotprod + &
                    & dipole_dist(frame2, i_at, i_coord) * &
                    & dipole_dist(frame2+frame1, i_at, i_coord)
          END DO

          sum_dip(layer) = sum_dip(layer) + second_leg(dotprod)

        ENDIF
      ENDDO
!$omp end parallel do

      DO i_bin = 1,nlayers
        IF ( counter_all(i_bin) /= 0 ) THEN
          surv_prob(i_bin) = surv_prob(i_bin) &
                           + dble(counter_in(i_bin)) / dble(counter_all(i_bin))
          !Averaged correlation function
          autocorr(frame1+1,i_bin) = autocorr(frame1+1,i_bin) + &
                                sum_dip(i_bin) / dble( counter_all(i_bin) ) 
        END IF
      END DO

    ENDDO

    surv_prob = surv_prob / dble(frame2-1) 

    DO i_bin = 1, nlayers

      IF ( surv_prob(i_bin) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(frame1+1,i_bin) = autocorr(frame1+1,i_bin) / dble( frame2-1 ) / surv_prob(i_bin)
      END IF

    END DO

  ENDDO

END SUBROUTINE ROT_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : autocorr
  USE parameters, ONLY : dt, iprint, box, nframes, nlayers
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION,PARAMETER :: au_to_fs = 2.418884326509D-2
  DOUBLE PRECISION           :: tc  

  tc = au_to_fs*dt

  OPEN(unit = 2,file = 'dipole_corr.dat')

  !Write the correlation function separated in bins
  write(2, fmt = '(a2,3x,*(E12.5,3x))'), "# ",&
            &( ( dble(i)-0.5 ) * box(3,3) / dble(nlayers), i=1,nlayers) 

  DO i = 0,nframes-1

    if (ANY(autocorr(1,:)==0.d0)) print *, "Something wrong"
    write(2,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*tc, &
         ( autocorr(i+1,j) / autocorr(1,j) , j=1,nlayers) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
