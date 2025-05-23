!Set of subroutines used to analyse the lipid trajectory

MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: oh_dist(:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
END MODULE histogram

PROGRAM OHcorr
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_POS
    CALL IDENTIFY_OH_GROUPS
    CALL OH_DISTRIBUTION (frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL ROT_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM OHcorr

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : autocorr, oh_dist
  USE parameters, ONLY : natoms, nframes, &
                         nOxygen, pos, coarse, &
                         nlayers, OH_index
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  CALL READ_INDEX ( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(coarse(nframes,nOxygen)); coarse=0
  ALLOCATE(oh_dist(3,2,nOxygen,nframes)); oh_dist = 0.d0
  ALLOCATE(autocorr(nlayers,nframes)); autocorr = 0
  ALLOCATE(OH_index(natoms)); OH_index =0

  OPEN(unit = 1,file = pos_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE OH_DISTRIBUTION (frame)
  !Build oh_dist matrix. oh_dist contains the normalized OH vectors 
  !for all OH groups in the system
  USE histogram,  ONLY : oh_dist 
  USE parameters, ONLY : natoms, coarse, nwater
  IMPLICIT NONE
  INTEGER                    :: frame, i, ih, layer
  INTEGER                    :: ind_H(2)
  INTEGER                    :: ind_layer
  INTEGER                    :: n_hydrogens, i_Ox, nOH
  INTEGER                    :: ihp 
  LOGICAL                    :: is_oxygen, switched_oh_index
  DOUBLE PRECISION           :: oh(2,3)

  !OH atom index
  nOH = 0
  i_Ox = 0 

  DO i = 1,natoms

    !Calculate the OH vector:
    IF ( is_oxygen(i) ) THEN

      i_Ox = i_Ox + 1
      CALL GET_IND_H(i, ind_H, n_hydrogens)

      IF ( n_hydrogens/=0 ) THEN 

        layer = ind_layer( i )

        CALL OH_VECT( i, ind_H, oh)

        coarse( frame, i_Ox ) = layer 

        DO ih = 1,n_hydrogens

          nOH = nOH + 1
          oh_dist(:,ih,i_Ox,frame) = oh(ih,:)
            
        END DO

      END IF !n_hydrogens

    END IF !is_oxygen

  ENDDO

  IF ( nOH /= 2*nwater ) THEN
    !print *, 'WARNING: Wrong number of OH groups ', nOH, ' on step ', frame
    !STOP
  END IF

END SUBROUTINE OH_DISTRIBUTION 

SUBROUTINE ROT_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : oh_dist, autocorr
  USE parameters, ONLY : nframes, nOxygen, coarse, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_dip(nlayers), surv_prob(nlayers)
  DOUBLE PRECISION           :: second_leg, dotprod 
  INTEGER                    :: frame1, frame2, iw, ipol
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers), ih 
  INTEGER                    :: layer
  INTEGER, PARAMETER         :: nsep=20
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,int(nframes/2.)

    surv_prob = 0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,nsep

      counter_in = 0
      counter_all = 0
      sum_dip = 0.D0

!$omp parallel do private(layer,dotprod) reduction(+:counter_in, sum_dip, counter_all)
      DO iw = 1,nOxygen

        layer = coarse(frame2,iw)
        IF (layer.eq.0) CYCLE
        counter_all(layer) = counter_all(layer) + 1

        !If the particle dont leave the window during the interval frame1
        IF ( in_layer(iw,frame2,frame2+frame1,50) ) THEN

          counter_in(layer) = counter_in(layer) + 1

          DO ih = 1,2

            dotprod = 0.d0

            DO ipol = 1,3

              dotprod = dotprod + &
                      & oh_dist(ipol, ih, iw, frame2) * &
                      & oh_dist(ipol, ih, iw, frame2+frame1)
            END DO

            sum_dip(layer) = sum_dip(layer) +  second_leg(dotprod)

          END DO

        ENDIF

      ENDDO
!$omp end parallel do

      DO layer = 1,nlayers
        IF ( counter_all(layer) /= 0 ) THEN
          surv_prob(layer) = surv_prob(layer) &
                           + dble(counter_in(layer)) / dble(counter_all(layer))
          !Averaged correlation function
          autocorr(layer,frame1+1) = autocorr(layer,frame1+1) + &
                                sum_dip(layer) / dble( counter_all(layer) ) 
        END IF
      END DO

    ENDDO

    surv_prob = surv_prob / dble(frame2-1) * dble(nsep)

    DO layer = 1, nlayers
      IF ( surv_prob(layer) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(layer,frame1+1) = autocorr(layer,frame1+1) / dble( frame2-1 ) / surv_prob(layer) * dble(nsep)
      END IF
    END DO

  ENDDO

END SUBROUTINE ROT_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : autocorr
  USE parameters, ONLY : nlayers, nframes, iprint, dt
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION           :: n_const(nlayers)
  DOUBLE PRECISION,PARAMETER :: au_to_fs = 2.418884326509D-2
  DOUBLE PRECISION           :: tc 

  tc = au_to_fs * dt

  OPEN( unit = 2, file = "OH_corr.dat" )

  DO i = 0,int(nframes/2.)

    write(2,fmt = '(E14.7,3X,*(E14.7,3X))'), dble(i)*tc, &
         ( autocorr(j,i+1) / autocorr(j,1) , j=1,nlayers) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
