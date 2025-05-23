!Set of subroutines used to analyse the lipid trajectory

MODULE histogram
  IMPLICIT NONE 
  REAL*8, ALLOCATABLE            :: oh_dist(:,:,:)
  REAL*8, ALLOCATABLE            :: oh_vel(:,:,:)
  REAL*8, ALLOCATABLE            :: autocorr(:,:)
  REAL*8, ALLOCATABLE            :: oh_pos(:,:,:)
  INTEGER                        :: tcorr, nOH
END MODULE histogram

PROGRAM SFG
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_EXTXYZ_POS_VEL
    CALL OH_DISTRIBUTION (frame)
  END DO

  CALL SSVAC 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM SFG

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : autocorr, oh_dist, oh_vel, &
                         oh_pos, tcorr, nOH
  USE parameters, ONLY : natoms, nframes, nwater, &
                         pos, vel, nequil, dt, atype, &
                         zoffset, nlayers, layers, coarse
  IMPLICIT NONE
  INTEGER                    :: i, get_number_atype
  CHARACTER(100)             :: pos_file, vel_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, dt, tcorr, nlayers
  READ(1,*) zoffset 
  IF(nlayers > 1) THEN
    ALLOCATE(layers(nlayers+1))
    READ(1, *) layers 
  END IF
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(vel(3,natoms)) ! velocities 
  ALLOCATE(atype(natoms))

  !Assuming all H are bound to O atoms
  OPEN(unit=1,file=pos_file)
  CALL READ_EXTXYZ_POS_VEL
  nOH = get_number_atype('H')
  CLOSE(1)

  ALLOCATE(oh_dist(3,nOH,nframes)); oh_dist = 0.d0
  ALLOCATE(oh_vel(3,nOH,nframes)); oh_vel = 0.d0
  ALLOCATE(oh_pos(3,nOH,nframes)); oh_pos = 0.d0
  ALLOCATE(coarse(nframes,nOH)); coarse = 1
  ALLOCATE(autocorr(nlayers,tcorr)); autocorr = 0

  OPEN(unit = 1,file = pos_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE OH_DISTRIBUTION (frame)
  !Compute OH_distance and OH velocity for each H atom
  USE histogram,  ONLY : oh_dist, oh_vel, oh_pos
  USE parameters, ONLY : natoms, pos, vel, atype, coarse
  IMPLICIT NONE
  INTEGER                    :: frame, i, ih
  INTEGER                    :: ind_OH, indO, layer_index
  REAL*8                     :: d_oh(3), sign_water_slab

  ind_OH=0
  DO i = 1,natoms
    IF ( trim(atype(i))=="H" ) THEN

      ind_OH = ind_OH + 1
      CALL MINDISTANCE_VECTOR_IND( i, "O", d_oh, indO ) 
      oh_dist (:,ind_OH,frame) = d_oh / norm2(d_oh) 
      oh_vel (:,ind_OH,frame) = vel(:,i) - vel(:,indO) 
      oh_pos (:,ind_OH,frame) = pos(:,indO) 
      coarse (frame,ind_OH) = layer_index(pos(3,indO),3)

    END IF

  ENDDO
 
END SUBROUTINE OH_DISTRIBUTION

SUBROUTINE SSVAC 
  !Surface specific velocity-velocity correlation function
  !Check eq. 9 of j. chem. phys 143, 124702 (2015)
  USE histogram,  ONLY : oh_dist, oh_vel, autocorr, &
                         tcorr, oh_pos, nOH
  USE parameters, ONLY : nframes, nwater, nlayers, coarse
  IMPLICIT NONE
  REAL*8                     :: surv_prob(nlayers)
  REAL*8                     :: alpha, alphamu, DistXYZ, d 
  REAL*8                     :: sign_water_slab
  REAL*8, PARAMETER          :: rcut=5.0d0
  INTEGER                    :: frame1, frame2, iw, jw, ipol, nsum
  INTEGER                    :: layer, i_bin, total_initial
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers) 
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    surv_prob = 0.d0
    nsum = 0
    counter_in = 0
    counter_all = 0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,100

      nsum = nsum + 1

      DO iw = 1, nOH

        layer = coarse(frame2,iw)

        IF ( layer .ne. 0 ) THEN

          counter_all(layer) = counter_all(layer) + 1

          IF ( in_layer(iw,frame2,frame2+frame1,1) ) THEN

            alpha = 0.d0
            alphamu = 0.d0
            counter_in(layer) = counter_in(layer) + 1

!$omp parallel do reduction(+:alpha) private(d)
            !DO jw = 1, nOH
            DO jw = iw, iw !Autocorrelation only
              d = DistXYZ(oh_pos(:,iw,frame2),oh_pos(:,jw,frame2))
              IF (d<rcut) THEN
                DO ipol = 1,3
                  alpha = alpha + &
                    & oh_vel(ipol, jw, frame2+frame1) * &
                    & oh_dist(ipol, jw, frame2+frame1)
                END DO ! ipol
              END IF ! d<rcut
            END DO ! jw
!$omp end parallel do
            alphamu = oh_vel(3, iw, frame2) * alpha !&
                                !* sign_water_slab(oh_pos(3,iw,frame2))

            autocorr(layer, frame1+1) = autocorr(layer, frame1+1) + alphamu 

          END IF ! in_layer
        END IF ! layer

      END DO ! iw
    ENDDO !frame2

    DO i_bin = 1,nlayers
      IF ( counter_all(i_bin) /= 0 ) THEN
        surv_prob(i_bin) =  &
                         dble(counter_in(i_bin)) / dble(counter_all(i_bin))
      END IF
    END DO

    !Averaged autocorrelation function
    DO i_bin = 1, nlayers
      IF ( surv_prob(i_bin) /= 0 ) THEN
        autocorr(i_bin,frame1+1) = autocorr(i_bin,frame1+1) / surv_prob(i_bin) / dble(nsum) 
      END IF
    END DO

  END DO ! frame1

END SUBROUTINE SSVAC 

SUBROUTINE PRINT_RESULTS 
  USE histogram,  ONLY : autocorr, tcorr
  USE parameters, ONLY : dt
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN( unit = 2, file = "ssvac_layers.dat" )

  DO i = 0,tcorr-1

    write(2,fmt = '(E13.6,*(3X,E14.7))') dble(i)*dt, autocorr(:,i+1) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

REAL*8 FUNCTION sign_water_slab(z)
  use parameters, only : box, zoffset
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: z
  REAL*8 :: posz

  posz = z + zoffset
  posz = posz - nint(posz/box(3,3))*box(3,3) 
  if (posz<0.0) then
      sign_water_slab = 1.d0
  else
      sign_water_slab = -1.d0
  end if

END FUNCTION sign_water_slab

