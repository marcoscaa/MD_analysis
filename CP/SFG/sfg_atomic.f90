!SFG trhough ssVAC. Taylor expanding polarizability and polarization
!as a function of the atomic coordinates. We always use the lab frame.
!Reading pre-computed derivatives from file "derivatives_atomic"

MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: h2o_vel(:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: der_pol(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: der_polbty(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: layerspf(:)
  INTEGER, ALLOCATABLE                     :: coarsepf(:,:) 
  INTEGER                                  :: nlayerspf
END MODULE histogram

PROGRAM SFG
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_POS_VEL
    CALL COARSE_GRAIN_POS (frame)
    CALL COARSE_GRAIN_PF (frame)
    CALL OH_DISTRIBUTION (frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL SSVAC 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM SFG

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : autocorr, h2o_vel, &
                         der_pol, der_polbty
  USE parameters, ONLY : natoms, nframes, nwater, &
                         pos, vel, coarse, nlayers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, vel_file, index_file,&
                                deriv_files

  CALL getarg(1, pos_file)
  CALL getarg(2, vel_file)
  CALL getarg(3, index_file)
  deriv_files="derivatives_atomic"

  CALL READ_INDEX ( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(vel(natoms,3)) ! velocities 
  ALLOCATE(coarse(nframes,nwater)); coarse=1
  ALLOCATE(h2o_vel(3,3,nwater,nframes)); h2o_vel = 0.d0
  ALLOCATE(autocorr(nlayers,nframes)); autocorr = 0

  CALL READ_DERIVATIVES ( deriv_files )

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = vel_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE OH_DISTRIBUTION (frame)
  !Generate the OH distribution for each water in a specific frame
  USE histogram,  ONLY : h2o_vel
  USE parameters, ONLY : natoms, vel
  IMPLICIT NONE
  INTEGER                    :: frame, i, ih
  INTEGER                    :: ind_O, ind_H(2)
  LOGICAL                    :: is_water_oxygen
  DOUBLE PRECISION           :: oh(3)

  !Water oxygen atom index
  ind_O = 0

  DO i = 1,natoms

    !Calculate the OH vector:
    IF ( is_water_oxygen(i) ) THEN

      ind_O = ind_O + 1
      ind_H = 0
      h2o_vel (:,1,ind_O,frame) = vel(i,:)
      CALL get_H(i, ind_H)
   
      !Two OH vectors, taking PBC into account 
      DO ih = 1,2

        h2o_vel (:,ih+1,ind_O,frame) = vel(ind_H(ih),:)

      END DO

    END IF

  ENDDO

END SUBROUTINE OH_DISTRIBUTION

SUBROUTINE SSVAC 
  !Surface specific velocity-velocity correlation function
  !Only the atomic auto-correlation function is considered for the moment
  USE histogram,  ONLY : h2o_vel, autocorr, &
                         der_pol, der_polbty, coarsepf
  USE parameters, ONLY : nframes, nwater, coarse, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_dip(nlayers), surv_prob(nlayers)
  DOUBLE PRECISION           :: dotprod(2) 
  INTEGER                    :: frame1, frame2, iw, ipol
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers)
  INTEGER                    :: ih, ih2
  INTEGER                    :: layer, layerpf
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,(nframes-1)/5.

    surv_prob = 0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,10

      counter_in = 0
      counter_all = 0
      sum_dip = 0.D0

!$omp parallel do private(layer,layerpf,dotprod) reduction(+:counter_in, sum_dip, counter_all)
      DO iw = 1,nwater

        layer = coarse(frame2,iw)
        layerpf = coarsepf(frame2,iw)
        counter_all(layer) = counter_all(layer) + 1

        !If the particle dont leave the window during the interval frame1
        IF ( in_layer(iw,frame2,frame2+frame1,10) ) THEN

          counter_in(layer) = counter_in(layer) + 1

          DO ih = 1,3
    
            DO ih2 = 1,3

              dotprod = 0.d0

              DO ipol = 1,3

                dotprod(1) = dotprod(1) + &
                             h2o_vel(ipol, ih, iw, frame2+frame1) * &
                             der_polbty(ipol,ih, layerpf) 

                dotprod(2) = dotprod(2) + &
                             h2o_vel(ipol, ih2, iw, frame2) * &
                             der_pol(ipol,ih2, layerpf)
              END DO

              dotprod(1) = dotprod(1) * dotprod(2)
              sum_dip(layer) = sum_dip(layer) + dotprod(1) 

            END DO

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

    surv_prob = surv_prob / dble(frame2-1) / 10.d0 

    DO layer = 1, nlayers
      IF ( surv_prob(layer) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(layer,frame1+1) = autocorr(layer,frame1+1) / dble( frame2-1 ) / surv_prob(layer) / 10.d0
      END IF
    END DO

  ENDDO

END SUBROUTINE SSVAC 

SUBROUTINE READ_DERIVATIVES ( derivatives_file )
  !Read pre-computed derivatives of polarization and polarizability
  !w.r.t. the rOH vector
  USE histogram,  ONLY : der_pol, der_polbty, layerspf, nlayerspf, &
                         coarsepf
  USE parameters, ONLY : nlayers, nframes, nwater
  CHARACTER(100), INTENT(IN) :: derivatives_file 
  INTEGER                    :: ipol, il, io
  INTEGER                    :: ilayer, itype
  DOUBLE PRECISION           :: tmp1(3), tmp2(3)

  OPEN(unit = 50,file = derivatives_file)

  READ(50,*), nlayerspf

  !We use the full layer definition here
  ALLOCATE(der_polbty(3,3,nlayerspf)); der_polbty = 0
  ALLOCATE(der_pol(3,3,nlayerspf)); der_pol = 0
  ALLOCATE(layerspf(nlayerspf+1)); layerspf = 0
  ALLOCATE(coarsepf(nframes,nwater)); coarsepf = 0

  READ(50,*), layerspf

  DO 

    READ(50,*,IOSTAT=io), ilayer, itype, tmp1, tmp2 

    IF ( io .eq. 0 ) THEN
      der_pol(:,itype,ilayer) = tmp1
      der_polbty(:,itype,ilayer) = tmp2
    ELSE
      EXIT
    END IF

  END DO
  
  !indexes 2 and 3 correspond to H atoms of the same water molecule
  der_pol(:,3,:) = der_pol(:,2,:) 
  der_polbty(:,3,:) = der_polbty(:,2,:)

  CLOSE(50)

END SUBROUTINE READ_DERIVATIVES

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : autocorr
  USE parameters, ONLY : nlayers, nframes, iprint, dt
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION           :: n_const(nlayers)
  DOUBLE PRECISION,PARAMETER :: au_to_fs = 2.418884326509D-2
  DOUBLE PRECISION,PARAMETER :: tc = au_to_fs*dt !*dble(iprint) 

  OPEN( unit = 2, file = "ssvac.dat" )

  DO i = 0,nframes-1

    write(2,fmt = '(E13.6,3X,*(E14.7,3X))'), dble(i)*tc, &
         ( autocorr(j,i+1) / autocorr(j,1), j=1,nlayers) 
         !( autocorr(j,i+1) , j=1,nlayers) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE COARSE_GRAIN_PF (frame)
  USE histogram, ONLY : nlayerspf, layerspf, coarsepf
  USE parameters, ONLY : natoms, box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, ind_z, iwat
  LOGICAL                              :: is_water_oxygen
  DOUBLE PRECISION                     :: z
  
  coarsepf(frame,:) = 0
  iwat = 0

  DO iat = 1, natoms

    !Selecting only OW
    IF ( is_water_oxygen(iat) ) THEN

      iwat = iwat + 1
      ind_z = 0
      z = pos(iat,3) - nint(pos(iat,3)/box(3,3))*box(3,3)
      z = z + box(3,3)/2.d0
      DO WHILE ( z > layerspf(ind_z) )
        ind_z = ind_z + 1
      END DO
      coarsepf(frame,iwat) = ind_z

    END IF

  END DO

END SUBROUTINE COARSE_GRAIN_PF
