!Set of subroutines used to analyse the lipid trajectory

MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: oh_dist(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: oh_vel(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: prefactor(:)
  DOUBLE PRECISION, ALLOCATABLE            :: layerspf(:)
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
    CALL OH_DISTRIBUTION (frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL SSVAC 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM SFG

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : autocorr, oh_dist, oh_vel,&
                         prefactor
  USE parameters, ONLY : natoms, nframes, nwater, &
                         pos, vel, coarse, nlayers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, vel_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, vel_file)
  CALL getarg(3, index_file)

  CALL READ_INDEX ( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(vel(natoms,3)) ! velocities 
  ALLOCATE(coarse(nframes,2*nwater)); coarse=0
  ALLOCATE(oh_dist(3,2*nwater,nframes)); oh_dist = 0.d0
  ALLOCATE(oh_vel(3,2*nwater,nframes)); oh_vel = 0.d0
  ALLOCATE(autocorr(nlayers,nframes)); autocorr = 0

  CALL READ_DERIVATIVES ( "derivatives_OH" )

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = vel_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE OH_DISTRIBUTION (frame)
  !Generate the OH distribution for each water in a specific frame
  USE histogram,  ONLY : oh_dist, oh_vel, prefactor
  USE parameters, ONLY : natoms, vel, coarse
  IMPLICIT NONE
  INTEGER                    :: frame, i, ih
  INTEGER                    :: ind_OH, ind_H(2)
  INTEGER                    :: ind_layer_pf
  INTEGER                    :: ind_layer
  LOGICAL                    :: is_oxygen
  DOUBLE PRECISION           :: oh(3)

  !OH bond index
  ind_OH = 0

  DO i = 1,natoms

    !Calculate the OH vector:
    IF ( is_oxygen(i) ) THEN

      ind_H = 0
      CALL get_H(i, ind_H)
   
      !Max two OH vectors per O atom, taking PBC into account 
      DO ih = 1,2

        IF ( ind_H(ih) .ne. 0 ) THEN 

          ind_OH = ind_OH + 1
          CALL OH_VECT( i, ind_H(ih), oh )
          oh_dist (:,ind_OH,frame) = oh(:) / norm2( oh(:) ) * &
                                      prefactor(ind_layer_pf(i))
          oh_vel (:,ind_OH,frame) = vel(ind_H(ih),:) - vel(i,:) 
          coarse(frame,ind_OH) = ind_layer(i)

        END IF

      END DO

    END IF

  ENDDO

END SUBROUTINE OH_DISTRIBUTION

SUBROUTINE SSVAC 
  !Surface specific velocity-velocity correlation function
  !Check eq. 9 of j. chem. phys 143, 124702 (2015)
  USE histogram,  ONLY : oh_dist, oh_vel, autocorr
  USE parameters, ONLY : nframes, nwater, coarse, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_dip(nlayers), surv_prob(nlayers)
  DOUBLE PRECISION           :: second_leg, dotprod 
  INTEGER                    :: frame1, frame2, iw, ipol
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers), ih 
  INTEGER                    :: layer
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,(nframes-1)/5.

    surv_prob = 0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,10

      counter_in = 0
      counter_all = 0
      sum_dip = 0.D0

!$omp parallel do private(layer,dotprod) reduction(+:counter_in, sum_dip, counter_all)
      DO iw = 1,2*nwater

        layer = coarse(frame2,iw)
        counter_all(layer) = counter_all(layer) + 1

        !If the particle dont leave the window during the interval frame1
        IF ( in_layer(iw,frame2,frame2+frame1,10) ) THEN

          counter_in(layer) = counter_in(layer) + 1

          dotprod = 0.d0

          DO ipol = 1,3

            dotprod = dotprod + &
                    & oh_vel(ipol, iw, frame2+frame1) * &
                    & oh_dist(ipol, iw, frame2+frame1)
          END DO

          dotprod = dotprod * oh_vel(3, iw, frame2)
          sum_dip(layer) = sum_dip(layer) + dotprod 

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
         ( autocorr(j,i+1) / autocorr(j,1) , j=1,nlayers) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE READ_DERIVATIVES( der_file )
  !User should provide file with the product of polarization 
  !and polarizability derivatives for each layer 
  !First line of the file: layers coordinates (z)
  USE histogram, ONLY : prefactor, nlayerspf, layerspf
  IMPLICIT NONE
  CHARACTER(14),INTENT(IN)   :: der_file

  OPEN(unit = 50,file = der_file)

  READ(50,*), nlayerspf

  ALLOCATE(prefactor(nlayerspf)); prefactor = 0
  ALLOCATE(layerspf(nlayerspf+1)); layerspf = 0

  READ(50,*), layerspf
  READ(50,*), prefactor
  
  CLOSE(50)

END SUBROUTINE READ_DERIVATIVES

SUBROUTINE OH_VECT(O, H, OH)
  ! Normalized vector from atom O to atom H
  USE parameters, ONLY : box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: O, H
  DOUBLE PRECISION, INTENT(INOUT)      :: OH(3)
  INTEGER                              :: ipol

  DO ipol=1,3

    OH(ipol) = pos(H,ipol) - pos(O,ipol)
    OH(ipol) = OH(ipol) - nint( OH(ipol)/box(ipol,ipol) ) * box(ipol,ipol)

  END DO

  OH = OH / norm2(OH)

END SUBROUTINE oh_vect

INTEGER FUNCTION ind_layer_pf( ind_O )
  !Layer index for prefactor. ind_O is the index of water oxygen in pos
  USE histogram, ONLY : nlayerspf, layerspf
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind_O
  INTEGER                    :: il
  DOUBLE PRECISION           :: z

  z = pos(ind_O,3)
  z = z - nint(z/box(3,3))*box(3,3)
  z = z + box(3,3)/2.d0

  DO il = 1,nlayerspf

    if ( z <= layerspf(il+1) ) THEN
      ind_layer_pf = il
      RETURN
    endif

  END DO

  ind_layer_pf = 0

END FUNCTION ind_layer_pf

INTEGER FUNCTION ind_layer( ind_O )
  !Layer index for prefactor. ind_O is the index of water oxygen in pos
  USE parameters, ONLY : pos, box, layers, nlayers
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind_O
  INTEGER                    :: il
  DOUBLE PRECISION           :: z

  z = pos(ind_O,3)
  z = z - nint(z/box(3,3))*box(3,3)
  z = z + box(3,3)/2.d0

  DO il = 1,nlayers

    if ( z <= layers(il+1) ) THEN
      ind_layer = il
      RETURN
    endif

  END DO

  ind_layer = 0

END FUNCTION ind_layer
